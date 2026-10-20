# Gene length check.
#
# The expected length of a gene is not given by the user: it is inferred
# from the gene sequences themselves, as the most frequent (modal) length
# among the OTHER unique sequences of the same gene in the same species
# (or, without a species column, the same comparison group). The focal
# sequence is left out so that it cannot define its own expectation.

"""
    modal_length(lengths) -> Union{Int,Nothing}

Most frequent value of `lengths`; `nothing` when empty or when two or more
lengths are equally frequent (no single expected length).
"""
function modal_length(lengths::AbstractVector{<:Integer})
    isempty(lengths) && return nothing
    c = Dict{Int,Int}()
    for x in lengths
        c[x] = get(c, x, 0) + 1
    end
    best = maximum(values(c))
    modes = [k for (k, v) in c if v == best]
    length(modes) == 1 ? only(modes) : nothing
end

"""
Where does the length difference of unique sequence `u` lie, relative to
the closest comparator of the expected length (at most 50 are examined)?

- `:terminal`: the sequences agree when anchored at the 5' end or at the
  3' end (mismatch fraction ≤ `threshold`), so the extra or missing bases
  are at an end;
- `:internal`: no end-anchored match, but the edit distance to the closest
  comparator is ≤ `threshold` × length, so the indel is inside the gene;
- `:undetermined`: no comparator close enough to decide.
"""
function _indel_location(g::GeneData, u::Int, comparators::Vector{Int}, threshold::Float64)
    isempty(comparators) && return :undetermined
    a = ungapped(g.masks[u])
    best_end = Inf; closest = 0
    for v in Iterators.take(comparators, 50)
        b = ungapped(g.masks[v])
        k = min(length(a), length(b))
        k == 0 && continue
        pre = hamming_pdistance(view(a, 1:k), view(b, 1:k))
        suf = hamming_pdistance(view(a, length(a) - k + 1:length(a)), view(b, length(b) - k + 1:length(b)))
        e = min(isnan(pre) ? Inf : pre, isnan(suf) ? Inf : suf)
        e < best_end && (best_end = e; closest = v)
    end
    best_end <= threshold && return :terminal
    closest == 0 && return :undetermined
    b = ungapped(g.masks[closest])
    L = max(length(a), length(b))
    kmax = floor(Int, threshold * L + 1e-9) + abs(length(a) - length(b))
    capped_edit_distance(a, b, kmax) <= kmax ? :internal : :undetermined
end

"""
Length test for every row and gene. Returns the `lengths` table and records
the `:length` test, `:length_outlier` (any difference larger than
`length_tolerance`) and `:invalid_frame` (difference not a multiple of 3).
"""
function length_analysis!(ctx::AnalysisContext, L::Ledger, allgroups::Dict{Symbol,Groups}, species_col)
    cfg = ctx.cfg; n = nrow(ctx.input)
    sp = _row_species(ctx, species_col)
    parts = DataFrame[]
    for g in ctx.genes
        gn = genename(g); G = allgroups[gn]
        scope(r) = sp === nothing ? G.row_cluster[r] : sp[r]
        scope_name = sp === nothing ? "comparison group" : "species"
        lens = [count(!=(GAP_MASK), m) for m in g.masks]
        members = Dict{Any,Vector{Int}}()
        for r in 1:n
            u = g.row_unique[r]; s = scope(r)
            (u == 0 || s === missing) && continue
            push!(get!(Vector{Int}, members, s), u)
        end
        for v in values(members)
            unique!(sort!(v))
        end
        I = Union{Missing,Int}
        len = Vector{I}(missing, n); modal = Vector{I}(missing, n); diff = Vector{I}(missing, n)
        ncomp = Vector{I}(missing, n); triplet = Vector{Union{Missing,Bool}}(missing, n)
        status = fill(:not_evaluated, n); reason = fill("", n)
        cache = Dict{Tuple{Any,Int},Any}()
        locs = Dict{Tuple{Any,Int},Symbol}()
        location = Vector{Union{Missing,Symbol}}(missing, n)
        for r in 1:n
            u = g.row_unique[r]
            if u == 0
                reason[r] = string(g.row_status[r]); continue
            end
            len[r] = lens[u]
            s = scope(r)
            if s === missing
                reason[r] = sp === nothing ? (G.row_reason[r] == "" ? "no comparison group" : G.row_reason[r]) :
                            "species missing: no within-species expected length"
                continue
            end
            m, k = get!(cache, (s, u)) do
                others = [lens[v] for v in members[s] if v != u]
                (length(others) >= cfg.length_min_comparators ? modal_length(others) : nothing, length(others))
            end
            ncomp[r] = k
            if k < cfg.length_min_comparators
                reason[r] = "$k other unique $gn sequence(s) in the $scope_name; at least " *
                            "$(cfg.length_min_comparators) needed to infer the expected length"
                continue
            elseif m === nothing
                reason[r] = "no single modal $gn length among the other sequences of the $scope_name"
                continue
            end
            d = lens[u] - m
            modal[r] = m; diff[r] = d; triplet[r] = d % 3 == 0; status[r] = :evaluated
            if abs(d) > cfg.length_tolerance
                loc = get!(locs, (s, u)) do
                    _indel_location(g, u, [v for v in members[s] if v != u && lens[v] == m],
                                    cfg.cluster_threshold)
                end
                location[r] = loc
                kind = d % 3 == 0 ? "a multiple of 3" : "not a multiple of 3"
                where_ = loc == :terminal ? "at a sequence end (terminal truncation or extension)" :
                         loc == :internal ? "inside the gene (internal indel)" :
                         "position undetermined (no closely related sequence of the expected length)"
                addflag!(L, r, gn, :length_outlier; metric = "length_difference", value = d,
                         threshold = cfg.length_tolerance,
                         detail = "$(lens[u]) nt vs expected $m nt (modal length of the other $k " *
                                  "unique $gn sequences of the $scope_name): difference $d nt, $kind, $where_")
                # Only an internal indel shifts the reading frame; a difference at
                # an end (e.g. an incomplete terminal codon) is not a frameshift.
                if d % 3 != 0 && loc == :internal
                    addflag!(L, r, gn, :invalid_frame; metric = "length_difference", value = d,
                             threshold = 0,
                             detail = "internal indel of $d nt relative to the closest sequence of the " *
                                      "expected $gn length: frameshift-like")
                end
            end
        end
        for r in 1:n
            status[r] == :evaluated ? addtest!(L, r, gn, :length, :evaluated) :
                                      addtest!(L, r, gn, :length, :not_evaluated; reason = reason[r])
        end
        push!(parts, DataFrame(row = 1:n, id = ctx.ids, gene = fill(gn, n),
                               scope = fill(sp === nothing ? :group : :species, n),
                               length_nt = len, modal_length_nt = modal, length_difference = diff,
                               multiple_of_3 = triplet, indel_location = location,
                               n_comparators = ncomp,
                               status = status, reason = reason))
    end
    reduce(vcat, parts)
end
