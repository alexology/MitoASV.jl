# ND2–COI (two-gene) concordance as a possible chimera indicator.
#
# Two genes of the same mitochondrial molecule share one genealogy, so ASVs
# that are close in one gene are expected to be close in the other. A PCR
# chimera joins a 5' part from one template to a 3' part from another, so
# its nearest neighbour can differ between the genes. The genes evolve at
# different rates, so only rank-based or neighbourhood comparisons are used
# and raw distances of the two genes are never compared directly.
# Discordance is a possible chimera indicator, not evidence: recombination
# is absent from animal mtDNA, but incomplete sampling, homoplasy, NUMTs and
# sequencing errors also create discordance.

"""
    spearman(x, y) -> Union{Missing,Float64}

Spearman rank correlation (average ranks for ties) over pairs where both
values are finite. `missing` with fewer than 3 pairs or zero rank variance.
"""
function spearman(x::AbstractVector, y::AbstractVector)
    keep = [i for i in eachindex(x, y) if _finite(x[i]) && _finite(y[i])]
    length(keep) < 3 && return missing
    rx = _ranks(Float64[x[i] for i in keep]); ry = _ranks(Float64[y[i] for i in keep])
    sx = std(rx); sy = std(ry)
    (sx == 0 || sy == 0) && return missing
    cor(rx, ry)
end
_finite(v) = v !== missing && isfinite(v)

function _ranks(v::Vector{Float64})
    n = length(v); p = sortperm(v); r = zeros(n)
    i = 1
    while i <= n
        j = i
        while j < n && v[p[j + 1]] == v[p[i]]
            j += 1
        end
        avg = (i + j) / 2
        for k in i:j
            r[p[k]] = avg
        end
        i = j + 1
    end
    r
end

"""Distance matrix (Float32, NaN = not comparable) among selected gene uniques."""
function _gene_matrix(items::Vector, dist; threaded::Bool)
    n = length(items)
    D = fill(NaN32, n, n)
    _tforeach(n; threaded) do i
        D[i, i] = 0f0
        for j in (i + 1):n
            d = dist(items[i], items[j])
            D[i, j] = isfinite(d) ? Float32(d) : NaN32
        end
    end
    for i in 1:n, j in (i + 1):n
        D[j, i] = D[i, j]
    end
    D
end

"""
    gene_concordance(seqs1, seqs2; k = 5, min_improvement = 0.005,
                     pair_threshold = 0.5, max_pairs = 1000, align_unequal = true,
                     proteins1 = nothing, proteins2 = nothing, abundance = nothing,
                     compare_trees = false, threaded = true) -> NamedTuple

Concordance between two genes of the same ASVs (`seqs1[i]` and `seqs2[i]`
come from ASV `i`; sequences as strings or mask vectors). Returns
`per_asv` (DataFrame), `pairs` (discordant pairs), `global_spearman_nt`,
`global_spearman_aa` and `tree_spearman`.

Per ASV:
- `nn_gene1`, `nn_gene2`: nearest neighbour in each gene (ties: lowest index);
- `knn_jaccard`: overlap of the `k` nearest neighbours in the two genes;
- `row_spearman_nt` (and `_aa`): rank correlation of the ASV's distances to
  all others in the two genes;
- `single_parent_distance`: best length-weighted distance to one other ASV
  over both genes; `two_parent_distance`: the same allowing a different
  ASV for each gene; `two_parent_improvement` = single − two (≥ 0). A
  large improvement with two different parents is the pattern expected for
  a chimera between them (as in UCHIME);
- `parents_abundance_ratio`: min(parent abundance) / ASV abundance, when
  abundances are given (chimeras are usually rarer than their parents).

Pairs are reported when the percentile ranks of their distances in the two
genes differ by ≥ `pair_threshold` (at most `max_pairs`, most discordant
first). Unequal-length sequences use a banded edit distance capped at
0.25; pairs beyond the cap are treated as not comparable.
"""
function gene_concordance(seqs1::AbstractVector, seqs2::AbstractVector; k::Integer = 5,
                          min_improvement::Real = 0.005, pair_threshold::Real = 0.5,
                          max_pairs::Integer = 1000, align_unequal::Bool = true,
                          proteins1 = nothing, proteins2 = nothing, abundance = nothing,
                          compare_trees::Bool = false, threaded::Bool = true)
    n = length(seqs1)
    n == length(seqs2) || throw(ArgumentError("both genes need one sequence per ASV"))
    m1 = [_maskof(s) for s in seqs1]; m2 = [_maskof(s) for s in seqs2]
    ntd(a, b) = sequence_distance(a, b; cutoff = length(a) == length(b) ? Inf : 0.25, align_unequal)
    D1 = _gene_matrix(m1, ntd; threaded); D2 = _gene_matrix(m2, ntd; threaded)
    A1 = proteins1 === nothing ? nothing : _gene_matrix(collect(proteins1), _protd; threaded)
    A2 = proteins2 === nothing ? nothing : _gene_matrix(collect(proteins2), _protd; threaded)
    l1 = [count(!=(GAP_MASK), m) for m in m1]; l2 = [count(!=(GAP_MASK), m) for m in m2]
    _concordance_core(D1, D2, A1, A2, l1, l2, abundance; k, min_improvement, pair_threshold,
                      max_pairs, compare_trees)
end

_maskof(s::Vector{UInt8}) = s
_maskof(s) = parse_sequence(s).mask
_protd(a, b) = (a === missing || b === missing) ? NaN : protein_pdistance(a, b)

function _concordance_core(D1, D2, A1, A2, l1, l2, abundance; k, min_improvement, pair_threshold,
                           max_pairs, compare_trees)
    n = size(D1, 1)
    F = Union{Missing,Float64}; I = Union{Missing,Int}
    nn1 = Vector{I}(missing, n); nn2 = Vector{I}(missing, n)
    d_nn1 = Vector{F}(missing, n); d_nn2 = Vector{F}(missing, n)
    jac = Vector{F}(missing, n); rs_nt = Vector{F}(missing, n); rs_aa = Vector{F}(missing, n)
    single = Vector{F}(missing, n); two = Vector{F}(missing, n); impr = Vector{F}(missing, n)
    ncand = zeros(Int, n); ratio = Vector{F}(missing, n)
    for i in 1:n
        cand = [j for j in 1:n if j != i && isfinite(D1[i, j]) && isfinite(D2[i, j])]
        ncand[i] = length(cand)
        isempty(cand) && continue
        w1 = l1[i] / (l1[i] + l2[i]); w2 = 1 - w1
        b1 = cand[argmin([(D1[i, j], j) for j in cand])]
        b2 = cand[argmin([(D2[i, j], j) for j in cand])]
        nn1[i] = b1; nn2[i] = b2
        d_nn1[i] = D1[i, b1]; d_nn2[i] = D2[i, b2]
        single[i] = minimum(w1 * D1[i, j] + w2 * D2[i, j] for j in cand)
        two[i] = w1 * D1[i, b1] + w2 * D2[i, b2]
        impr[i] = max(single[i] - two[i], 0.0)
        s1 = _knn_set(D1, i, cand, k); s2 = _knn_set(D2, i, cand, k)
        jac[i] = length(intersect(s1, s2)) / length(union(s1, s2))
        if length(cand) >= 5
            rs_nt[i] = spearman([D1[i, j] for j in cand], [D2[i, j] for j in cand])
            if A1 !== nothing && A2 !== nothing
                rs_aa[i] = spearman([A1[i, j] for j in cand], [A2[i, j] for j in cand])
            end
        end
        if abundance !== nothing && b1 != b2
            a = abundance[i]; p = min(coalesce(abundance[b1], NaN), coalesce(abundance[b2], NaN))
            (a !== missing && a > 0 && isfinite(p)) && (ratio[i] = p / a)
        end
    end
    flag = [impr[i] !== missing && impr[i] >= min_improvement && nn1[i] != nn2[i] for i in 1:n]
    # global and pairwise statistics
    iu = [(i, j) for i in 1:n for j in (i + 1):n]
    v1 = Float64[D1[i, j] for (i, j) in iu]; v2 = Float64[D2[i, j] for (i, j) in iu]
    gs_nt = spearman(v1, v2)
    gs_aa = (A1 === nothing || A2 === nothing) ? missing :
            spearman(Float64[A1[i, j] for (i, j) in iu], Float64[A2[i, j] for (i, j) in iu])
    ok = [isfinite(a) && isfinite(b) for (a, b) in zip(v1, v2)]
    pr = DataFrame(i = Int[], j = Int[], distance_gene1 = Float64[], distance_gene2 = Float64[],
                   percentile_gene1 = Float64[], percentile_gene2 = Float64[], discordance = Float64[])
    if count(ok) >= 3
        idx = findall(ok)
        p1 = _ranks(v1[idx]) ./ length(idx); p2 = _ranks(v2[idx]) ./ length(idx)
        disc = abs.(p1 .- p2)
        sel = findall(>=(pair_threshold), disc)
        sel = sel[sortperm([(-disc[s], idx[s]) for s in sel])]
        for s in sel[1:min(end, max_pairs)]
            (i, j) = iu[idx[s]]
            push!(pr, (i, j, v1[idx[s]], v2[idx[s]], p1[s], p2[s], disc[s]))
        end
    end
    tree = missing; tree_reason = "not requested"
    if compare_trees
        if all(isfinite, D1) && all(isfinite, D2) && n >= 3
            C1 = _cophenetic(Float64.(D1)); C2 = _cophenetic(Float64.(D2))
            tree = spearman([C1[i, j] for (i, j) in iu], [C2[i, j] for (i, j) in iu])
            tree_reason = "UPGMA trees per gene; Spearman correlation of cophenetic distances"
        else
            tree_reason = "not computed: some distances are not comparable or fewer than 3 ASVs"
        end
    end
    per = DataFrame(nn_gene1 = nn1, nn_gene2 = nn2, distance_nn_gene1 = d_nn1,
                    distance_nn_gene2 = d_nn2, knn_jaccard = jac, row_spearman_nt = rs_nt,
                    row_spearman_aa = rs_aa, single_parent_distance = single,
                    two_parent_distance = two, two_parent_improvement = impr,
                    parents_abundance_ratio = ratio, n_comparators = ncand,
                    discordant = flag)
    (per_asv = per, pairs = pr, global_spearman_nt = gs_nt, global_spearman_aa = gs_aa,
     tree_spearman = tree, tree_reason = tree_reason)
end

"""Cophenetic distances of the UPGMA (average linkage) tree."""
function _cophenetic(D::Matrix{Float64})
    n = size(D, 1)
    merges = nn_chain_merges(copy(D), :average)
    members = Dict(i => [i] for i in 1:n)
    C = zeros(n, n)
    # NN-chain does not emit merges in height order; sort them (heights are
    # monotone for average linkage) and replay with union-find.
    u = UnionFind(n)
    for (a, b, h) in sort(merges; by = x -> x[3])
        ra, rb = find!(u, a), find!(u, b)
        for x in members[ra], y in members[rb]
            C[x, y] = h; C[y, x] = h
        end
        uf_union!(u, ra, rb)
        r = find!(u, ra)
        members[r] = vcat(members[ra], members[rb])
        r == ra || delete!(members, ra)
        r == rb || delete!(members, rb)
    end
    C
end

# --- pipeline -------------------------------------------------------------------

const _CONC_FLOAT_COLS = (:distance_nn_gene1, :distance_nn_gene2, :knn_jaccard, :row_spearman_nt,
                          :row_spearman_aa, :single_parent_distance, :two_parent_distance,
                          :two_parent_improvement, :parents_abundance_ratio,
                          :global_spearman_nt, :global_spearman_aa, :tree_spearman)

"""
Two-gene concordance. With `concordance_scope = :species` (the default) and
a species column, ASVs are compared only with ASVs of the same species;
with `:all`, or without a species column, with every ASV of the dataset.
"""
function concordance_analysis!(ctx::AnalysisContext, L::Ledger, species_col = nothing)
    cfg = ctx.cfg; n = nrow(ctx.input)
    empty_pairs = DataFrame()
    sp = _row_species(ctx, species_col)
    notrows(status, reason) = DataFrame(row = 1:n, id = ctx.ids,
                                        status = fill(status, n), reason = fill(reason, n))
    if length(ctx.genes) < 2
        for r in 1:n
            addtest!(L, r, missing, :gene_concordance, :not_applicable; reason = "single-gene analysis")
        end
        return (table = notrows(:not_applicable, "single-gene analysis"), pairs = empty_pairs)
    end
    names = cfg.concordance_genes === nothing ? (genename(ctx.genes[1]), genename(ctx.genes[2])) :
            cfg.concordance_genes
    gi = [findfirst(g -> genename(g) == nm, ctx.genes) for nm in names]
    any(isnothing, gi) && throw(ArgumentError("concordance_genes $(names) are not genes of this analysis"))
    g1, g2 = ctx.genes[gi[1]], ctx.genes[gi[2]]
    tname = Symbol(names[1], "_", names[2])
    if !cfg.concordance
        for r in 1:n
            addtest!(L, r, missing, :gene_concordance, :not_evaluated; reason = "disabled (concordance = false)")
        end
        return (table = notrows(:not_evaluated, "disabled (concordance = false)"), pairs = empty_pairs)
    end
    by_species = sp !== nothing && cfg.concordance_scope == :species
    scope_of = [by_species ? sp[r] : 1 for r in 1:n]
    F = Union{Missing,Float64}
    cols = Dict{Symbol,Vector{F}}(c => Vector{F}(missing, n) for c in _CONC_FLOAT_COLS)
    nn1id = Vector{Any}(missing, n); nn2id = Vector{Any}(missing, n)
    ncomp = Vector{Union{Missing,Int}}(missing, n); disc = Vector{Union{Missing,Bool}}(missing, n)
    status = fill(:not_evaluated, n); reason = fill("", n)
    pairs_parts = DataFrame[]
    for r in 1:n
        if by_species && sp[r] === missing
            reason[r] = "species missing: no within-species comparison"
        elseif g1.row_unique[r] == 0 || g2.row_unique[r] == 0
            reason[r] = "both genes are needed"
        end
    end
    scopes = unique(scope_of[r] for r in 1:n if reason[r] == "")
    sname(s) = by_species ? "species $(s)" : "the dataset"
    for s in scopes
        rows = [r for r in 1:n if reason[r] == "" && isequal(scope_of[r], s)]
        unit_of = Dict{Int,Int}(); idx = Dict{Tuple{Int,Int},Int}(); units = Tuple{Int,Int}[]
        for r in rows
            key = (g1.row_unique[r], g2.row_unique[r])
            unit_of[r] = get!(idx, key) do
                push!(units, key); length(units)
            end
        end
        nu = length(units)
        if nu < 3 || nu > cfg.max_concordance_asvs
            why = nu < 3 ? "fewer than 3 distinct ASVs with both $(names[1]) and $(names[2]) in $(sname(s))" :
                  "$nu distinct ASVs in $(sname(s)) exceed max_concordance_asvs = $(cfg.max_concordance_asvs)"
            for r in rows
                reason[r] = why
            end
            continue
        end
        s1 = sort!(unique(first.(units))); s2 = sort!(unique(last.(units)))
        p1 = Dict(u => k for (k, u) in enumerate(s1)); p2 = Dict(u => k for (k, u) in enumerate(s2))
        ntd(a, b) = sequence_distance(a, b; cutoff = length(a) == length(b) ? Inf : 0.25,
                                      align_unequal = cfg.align_unequal)
        G1 = _gene_matrix([g1.masks[u] for u in s1], ntd; threaded = cfg.threaded)
        G2 = _gene_matrix([g2.masks[u] for u in s2], ntd; threaded = cfg.threaded)
        P1 = _gene_matrix([g1.proteins[u] for u in s1], _protd; threaded = cfg.threaded)
        P2 = _gene_matrix([g2.proteins[u] for u in s2], _protd; threaded = cfg.threaded)
        a = [p1[u[1]] for u in units]; b = [p2[u[2]] for u in units]
        l1 = [count(!=(GAP_MASK), g1.masks[u[1]]) for u in units]
        l2 = [count(!=(GAP_MASK), g2.masks[u[2]]) for u in units]
        ab = nothing
        if ctx.abundance !== nothing
            ab = Vector{Union{Missing,Float64}}(fill(0.0, nu))
            for r in rows
                x = ctx.abundance[r]; u = unit_of[r]
                ab[u] = (x === missing || ab[u] === missing) ? missing : ab[u] + x
            end
        end
        res = _concordance_core(G1[a, a], G2[b, b], P1[a, a], P2[b, b], l1, l2, ab;
                                k = cfg.concordance_k, min_improvement = cfg.chimera_min_improvement,
                                pair_threshold = cfg.pair_discordance_threshold,
                                max_pairs = cfg.max_discordant_pairs, compare_trees = cfg.compare_gene_trees)
        first_row = Dict{Int,Int}()
        for r in rows
            haskey(first_row, unit_of[r]) || (first_row[unit_of[r]] = r)
        end
        per = res.per_asv
        for r in rows
            u = unit_of[r]
            for c in _CONC_FLOAT_COLS[1:9]
                cols[c][r] = per[u, c]
            end
            cols[:global_spearman_nt][r] = res.global_spearman_nt
            cols[:global_spearman_aa][r] = res.global_spearman_aa
            cols[:tree_spearman][r] = res.tree_spearman
            x = per.nn_gene1[u]; nn1id[r] = x === missing ? missing : ctx.ids[first_row[x]]
            x = per.nn_gene2[u]; nn2id[r] = x === missing ? missing : ctx.ids[first_row[x]]
            ncomp[r] = per.n_comparators[u]; disc[r] = per.discordant[u]
            if per.n_comparators[u] > 0
                status[r] = :evaluated
            else
                reason[r] = "no comparable ASV in both genes in $(sname(s))"
            end
        end
        pr = res.pairs
        if nrow(pr) > 0
            push!(pairs_parts, DataFrame(
                species = Any[by_species ? s : missing for _ in 1:nrow(pr)],
                row_i = [first_row[i] for i in pr.i], id_i = [ctx.ids[first_row[i]] for i in pr.i],
                row_j = [first_row[j] for j in pr.j], id_j = [ctx.ids[first_row[j]] for j in pr.j],
                gene_pair = fill(tname, nrow(pr)),
                distance_gene1 = pr.distance_gene1, distance_gene2 = pr.distance_gene2,
                percentile_gene1 = pr.percentile_gene1, percentile_gene2 = pr.percentile_gene2,
                discordance = pr.discordance))
        end
    end
    tab = DataFrame(row = 1:n, id = ctx.ids, gene_pair = fill(tname, n),
                    scope = fill(by_species ? :species : :all, n),
                    species = sp === nothing ? fill(missing, n) : sp,
                    nn_gene1_id = nn1id, nn_gene2_id = nn2id)
    for c in _CONC_FLOAT_COLS[1:9]
        tab[!, c] = cols[c]
    end
    tab.n_comparators = ncomp; tab.discordant = disc
    for c in _CONC_FLOAT_COLS[10:12]
        tab[!, c] = cols[c]
    end
    tab.status = status; tab.reason = reason
    for r in 1:n
        if status[r] == :evaluated
            addtest!(L, r, missing, :gene_concordance, :evaluated)
            if disc[r] === true
                addflag!(L, r, missing, :gene_discordance; metric = "two_parent_improvement",
                         value = cols[:two_parent_improvement][r], threshold = cfg.chimera_min_improvement,
                         detail = "nearest $(names[1]) neighbour $(nn1id[r]) differs from nearest " *
                                  "$(names[2]) neighbour $(nn2id[r]); two parents fit better than one " *
                                  "by $(_r(cols[:two_parent_improvement][r])) (possible chimera indicator; " *
                                  "compared within $(by_species ? "the species" : "the whole dataset"))")
            end
        else
            addtest!(L, r, missing, :gene_concordance, :not_evaluated; reason = reason[r])
        end
    end
    (table = tab, pairs = isempty(pairs_parts) ? empty_pairs : reduce(vcat, pairs_parts))
end

"""The `k` nearest candidates of `i`, plus any candidate tied with the k-th distance."""
function _knn_set(D, i, cand, k)
    d = sort([D[i, j] for j in cand])
    dk = d[min(k, length(d))]
    Set(j for j in cand if D[i, j] <= dk)
end
