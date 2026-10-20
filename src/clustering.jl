# Reference-free grouping of ASVs into comparison groups.
#
# Connected components link two sequences whenever their distance is ≤ the
# threshold. Components can chain: A–B and B–C at 3% put A and C in one
# group even if A and C differ by 6%. Hierarchical clustering with average
# or complete linkage limits chaining but needs a full distance matrix.
# The default 3% threshold is a configurable example, not a species boundary.

mutable struct UnionFind
    parent::Vector{Int}
end
UnionFind(n::Int) = UnionFind(collect(1:n))
function find!(u::UnionFind, x::Int)
    while u.parent[x] != x
        u.parent[x] = u.parent[u.parent[x]]
        x = u.parent[x]
    end
    x
end
function uf_union!(u::UnionFind, a::Int, b::Int)
    ra, rb = find!(u, a), find!(u, b)
    ra == rb && return
    ra < rb ? (u.parent[rb] = ra) : (u.parent[ra] = rb)
end

"""Cluster labels 1..K numbered by the first member in input order."""
function _labels(u::UnionFind, n::Int)
    lab = zeros(Int, n); seen = Dict{Int,Int}()
    for i in 1:n
        r = find!(u, i)
        lab[i] = get!(seen, r, length(seen) + 1)
    end
    lab
end

"""
Edges (i < j) with distance ≤ threshold, computed in parallel. Each task
returns its own buffer; buffers are merged and sorted, so the output does
not depend on scheduling. No distance matrix is stored.
"""
function threshold_edges(dist, n::Int, threshold::Float64; threaded::Bool = true)
    nchunks = (threaded && Threads.nthreads() > 1) ? min(n, 8 * Threads.nthreads()) : 1
    work(c) = begin
        buf = Tuple{Int,Int,Float64}[]
        for i in c:nchunks:n, j in (i + 1):n
            d = dist(i, j, threshold)
            (!isnan(d) && d <= threshold) && push!(buf, (i, j, d))
        end
        buf
    end
    bufs = nchunks == 1 ? [work(1)] : fetch.([Threads.@spawn(work(c)) for c in 1:nchunks])
    edges = reduce(vcat, bufs; init = Tuple{Int,Int,Float64}[])
    sort!(edges)
end

function _components(dist, n, threshold; threaded = true)
    edges = threshold_edges(dist, n, threshold; threaded)
    u = UnionFind(n)
    nn = fill(Inf, n)
    for (i, j, d) in edges
        uf_union!(u, i, j)
        nn[i] = min(nn[i], d); nn[j] = min(nn[j], d)
    end
    _labels(u, n), nn
end

"""
Full distance matrix (diagonal 0). Pairs that are not comparable get `Inf`
so that they are never merged.
"""
function distance_matrix(dist, n::Int; cutoff::Float64 = Inf, threaded::Bool = true)
    D = zeros(n, n)
    _tforeach(n; threaded) do i
        for j in (i + 1):n
            d = dist(i, j, cutoff)
            D[i, j] = isnan(d) ? Inf : (isinf(d) ? max(cutoff, 0.0) : d)
        end
    end
    for i in 1:n, j in (i + 1):n
        D[j, i] = D[i, j]
    end
    D
end

"""
Agglomerative clustering with the nearest-neighbour chain algorithm
(Lance-Williams updates; `:single`, `:complete`, `:average`). Returns the
merges as `(a, b, height)` with `a`, `b` original indices representing
each cluster. `D` is overwritten.
"""
function nn_chain_merges(D::Matrix{Float64}, linkage::Symbol)
    n = size(D, 1)
    active = trues(n); sz = ones(Int, n)
    merges = Tuple{Int,Int,Float64}[]
    chain = Int[]
    nactive = n
    while nactive > 1
        isempty(chain) && push!(chain, findfirst(active))
        a = chain[end]
        prevc = length(chain) >= 2 ? chain[end - 1] : 0
        b = 0; best = Inf
        if prevc != 0
            b = prevc; best = D[a, prevc]
        end
        for j in 1:n
            (active[j] && j != a) || continue
            if D[a, j] < best || b == 0
                best = D[a, j]; b = j
            end
        end
        if b == prevc
            pop!(chain); pop!(chain)
            lo, hi = minmax(a, b)
            push!(merges, (lo, hi, best))
            for k in 1:n
                (active[k] && k != lo && k != hi) || continue
                dl, dh = D[lo, k], D[hi, k]
                v = linkage == :single ? min(dl, dh) :
                    linkage == :complete ? max(dl, dh) :
                    (sz[lo] * dl + sz[hi] * dh) / (sz[lo] + sz[hi])
                D[lo, k] = v; D[k, lo] = v
            end
            sz[lo] += sz[hi]; active[hi] = false; nactive -= 1
            filter!(x -> x != hi && x != lo, chain)
        else
            push!(chain, b)
        end
    end
    merges
end

function _hierarchical(dist, n, threshold, linkage; threaded = true)
    cap = max(0.5, 2threshold)
    D = distance_matrix(dist, n; cutoff = cap, threaded)
    nn = [n == 1 ? Inf : minimum(D[i, j] for j in 1:n if j != i) for i in 1:n]
    nn = [x <= threshold ? x : Inf for x in nn]
    merges = nn_chain_merges(D, linkage)
    u = UnionFind(n)
    for (a, b, h) in merges
        h <= threshold && uf_union!(u, a, b)
    end
    _labels(u, n), nn
end

"""
    cluster_asvs(sequences; method = :connected_components, threshold = 0.03,
                 level = :nucleotide, linkage = :average, align_unequal = true,
                 genetic_code = 5, frame = 1, threaded = true) -> DataFrame

Group sequences by pairwise distance without any reference.

- `method = :connected_components`: single-linkage components of the graph
  with an edge for every distance ≤ `threshold`. Only edges are stored
  (no distance matrix). Components can chain relatively divergent
  sequences through intermediates.
- `method = :hierarchical`: agglomerative clustering (`linkage =
  :average`, `:complete` or `:single`) cut at `threshold`; needs an n × n
  matrix. Distances above `max(0.5, 2 × threshold)` are stored at that cap.
- `level = :nucleotide`: p-distance on equal-length sequences (pairwise
  deletion of gaps and ambiguous bases); unequal lengths use a banded edit
  distance (`align_unequal = true`) or are not linked.
- `level = :protein`: p-distance between translations (`genetic_code`,
  `frame`), equal lengths only.

Identical sequences always share a cluster. `missing` sequences get a
`missing` cluster. The 3% threshold is an example, not a universal
species boundary. Returns `index`, `cluster`, `cluster_size` (number of
input sequences in the cluster) and `nn_distance` (distance to the nearest
other unique sequence if ≤ threshold, else `missing`).
"""
function cluster_asvs(sequences::AbstractVector; method::Symbol = :connected_components,
                      threshold::Real = 0.03, level::Symbol = :nucleotide,
                      linkage::Symbol = :average, align_unequal::Bool = true,
                      genetic_code = 5, frame::Integer = 1, threaded::Bool = true)
    method in (:connected_components, :hierarchical) ||
        throw(ArgumentError("method must be :connected_components or :hierarchical"))
    level in (:nucleotide, :protein) || throw(ArgumentError("level must be :nucleotide or :protein"))
    n = length(sequences)
    unit = zeros(Int, n); index = Dict{Any,Int}(); items = Any[]
    code = MitoASV.genetic_code(genetic_code)
    for i in 1:n
        p = parse_sequence(sequences[i])
        (p === missing || isempty(p.mask)) && continue
        key = level == :nucleotide ? p.mask : read_frame(p.mask, Int(frame), code).protein
        unit[i] = get!(index, key) do
            push!(items, key); length(items)
        end
    end
    labels, nn = cluster_units(items; method, threshold = Float64(threshold), level, linkage,
                               align_unequal, threaded)
    cl = Union{Missing,Int}[unit[i] == 0 ? missing : labels[unit[i]] for i in 1:n]
    sizes = Dict{Int,Int}()
    for c in skipmissing(cl)
        sizes[c] = get(sizes, c, 0) + 1
    end
    DataFrame(index = 1:n, cluster = cl,
              cluster_size = Union{Missing,Int}[c === missing ? missing : sizes[c] for c in cl],
              nn_distance = Union{Missing,Float64}[unit[i] == 0 || !isfinite(nn[unit[i]]) ? missing :
                                                   nn[unit[i]] for i in 1:n])
end

"""
Cluster unique items (mask vectors or protein strings). Returns labels and
nearest-neighbour distances (`Inf` when no neighbour is within threshold).
"""
function cluster_units(items::AbstractVector; method::Symbol, threshold::Float64, level::Symbol,
                       linkage::Symbol = :average, align_unequal::Bool = true,
                       threaded::Bool = true, max_hierarchical::Int = typemax(Int))
    n = length(items)
    n == 0 && return Int[], Float64[]
    dist = if level == :nucleotide
        (i, j, c) -> sequence_distance(items[i], items[j]; cutoff = c, align_unequal)
    else
        (i, j, c) -> protein_pdistance(items[i], items[j]; cutoff = c)
    end
    if method == :hierarchical
        n > max_hierarchical && throw(ArgumentError(
            "hierarchical clustering of $n unique sequences exceeds max_hierarchical = " *
            "$max_hierarchical (it needs an n × n matrix); use :connected_components"))
        return _hierarchical(dist, n, threshold, linkage; threaded)
    end
    _components(dist, n, threshold; threaded)
end
