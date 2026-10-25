# Comparison groups and within-group comparative analyses.

"""Sequence used to group each row (`nothing` when unavailable)."""
function _cluster_keys(ctx::AnalysisContext)
    cfg = ctx.cfg; n = nrow(ctx.input)
    keys = Vector{Any}(nothing, n)
    if cfg.cluster_on !== nothing
        k = findfirst(g -> genename(g) == cfg.cluster_on, ctx.genes)
        k === nothing && throw(ArgumentError("cluster_on = $(cfg.cluster_on) is not a gene of this analysis"))
        keys, _ = _gene_keys(ctx, ctx.genes[k])
        basis = cfg.cluster_level == :protein ? "protein translation" : "nucleotide sequence"
        return keys, "$basis of gene $(cfg.cluster_on)"
    end
    if cfg.cluster_level == :protein
        for r in 1:n
            parts = String[]
            for g in ctx.genes
                u = g.row_unique[r]
                (u == 0 || g.proteins[u] === missing) && (parts = String[]; break)
                push!(parts, g.proteins[u])
            end
            isempty(parts) || (keys[r] = join(parts, "|"))
        end
        return keys, "protein translation of " * join(genename.(ctx.genes), "+")
    end
    if ctx.main !== nothing
        for r in 1:n
            p = ctx.main[r]
            (p === missing || isempty(p.mask)) || (keys[r] = p.mask)
        end
        return keys, "full input sequence"
    end
    for r in 1:n
        if all(g -> g.row_unique[r] != 0, ctx.genes)
            keys[r] = reduce(vcat, [g.masks[g.row_unique[r]] for g in ctx.genes])
        end
    end
    keys, "concatenated gene sequences"
end

struct Groups
    row_cluster::Vector{Union{Missing,Int}}
    labels::Vector{Any}
    group_species::Vector{Any}          # species of each group (`missing` without species_col)
    row_reason::Vector{String}          # why a row has no group ("" if it has one)
    source::Symbol
    table::DataFrame
end

"""
Species of each row, or `nothing` without `species_col`. Comparison groups
never contain more than one species.
"""
_row_species(ctx::AnalysisContext, species_col) =
    species_col === nothing ? nothing : collect(ctx.input[!, species_col])

"""Grouping sequence of each row for one gene (`nothing` when unavailable)."""
function _gene_keys(ctx::AnalysisContext, g::GeneData)
    n = nrow(ctx.input)
    keys = Vector{Any}(nothing, n)
    for r in 1:n
        u = g.row_unique[r]
        u == 0 && continue
        if ctx.cfg.cluster_level == :protein
            g.proteins[u] === missing || (keys[r] = g.proteins[u])
        else
            keys[r] = g.masks[u]
        end
    end
    keys, (ctx.cfg.cluster_level == :protein ? "protein" : "nucleotide sequence") *
          " of gene $(genename(g)) (grouped separately for each gene)"
end

"""
Comparison groups for every gene, as a `Dict` gene => `Groups`.

- With a full-sequence column (`sequence_col`) and `cluster_on = nothing`, or
  with `cluster_on = <gene>`, all genes share the same groups.
- Without a full-sequence column (one column per gene), or with
  `cluster_on = :per_gene`, each gene is grouped on its own sequences, so a
  row missing one gene still has groups for the others.
Groups never mix species when `species_col` is given.
"""
function assign_groups!(ctx::AnalysisContext, L::Ledger, cluster_col, species_col = nothing)
    cfg = ctx.cfg
    per_gene = cfg.cluster_on === :per_gene || (cfg.cluster_on === nothing && ctx.main === nothing)
    out = Dict{Symbol,Groups}()
    if per_gene
        for (k, g) in enumerate(ctx.genes)
            keys, basis = _gene_keys(ctx, g)
            out[genename(g)] = _assign_groups(ctx, L, keys, basis, cluster_col, species_col, k == 1)
        end
    else
        keys, basis = _cluster_keys(ctx)
        G = _assign_groups(ctx, L, keys, basis, cluster_col, species_col, true)
        for g in ctx.genes
            out[genename(g)] = G
        end
    end
    out
end

"""The `clusters` table: one row per input row and gene."""
function clusters_table(groups::Dict{Symbol,Groups}, gds)
    reduce(vcat, [insertcols!(copy(groups[genename(g)].table), 3, :gene => fill(genename(g), nrow(groups[genename(g)].table)))
                  for g in gds])
end

function _assign_groups(ctx::AnalysisContext, L::Ledger, keys, basis, cluster_col, species_col,
                        flag_species::Bool)
    cfg = ctx.cfg; n = nrow(ctx.input)
    sp = _row_species(ctx, species_col)
    scope(r) = sp === nothing ? 1 : sp[r]
    row_cluster = Vector{Union{Missing,Int}}(missing, n)
    nn = Vector{Union{Missing,Float64}}(missing, n)
    unit = Vector{Union{Missing,Int}}(missing, n)
    row_reason = fill("", n)
    labels = Any[]; gspecies = Any[]
    # Units are distinct (species, sequence) pairs: a sequence assigned to two
    # species is analysed separately in each.
    ukeys = Vector{Any}(nothing, n)
    for r in 1:n
        if sp !== nothing && sp[r] === missing
            row_reason[r] = "species missing: no within-species comparison group"
            flag_species && addflag!(L, r, missing, :missing_species;
                                     detail = "value of $(species_col) is missing")
            continue
        end
        keys[r] === nothing && (row_reason[r] = "no sequence to group"; continue)
        ukeys[r] = (scope(r), keys[r])
    end
    items = _units!(unit, ukeys)
    unit_scope = [first(it) for it in items]
    scopes = unique(unit_scope)
    newgroup!(label, s) = (push!(labels, label); push!(gspecies, sp === nothing ? missing : s); length(labels))
    if cluster_col !== nothing
        source = :user
        basis = "user column $(cluster_col)"
        idx = Dict{Any,Int}()
        for (r, v) in enumerate(ctx.input[!, cluster_col])
            unit[r] === missing && continue
            if v === missing
                row_reason[r] = "user group ($(cluster_col)) missing"; continue
            end
            key = sp === nothing ? v : (sp[r], v)
            row_cluster[r] = get!(() -> newgroup!(key, scope(r)), idx, key)
        end
    elseif !cfg.clustering || cfg.cluster_method == :none
        source = :single_group
        basis *= sp === nothing ? " (clustering disabled: all sequences form one comparison group)" :
                 " (clustering disabled: each species forms one comparison group)"
        idx = Dict{Any,Int}()
        for r in 1:n
            unit[r] === missing && continue
            row_cluster[r] = get!(() -> newgroup!(scope(r), scope(r)), idx, scope(r))
        end
    else
        source = cfg.cluster_method
        sp === nothing || (basis *= ", clustered within each species")
        gid = Dict{Tuple{Any,Int},Int}()
        uclust = zeros(Int, length(items)); unn = fill(Inf, length(items))
        for s in scopes
            us = findall(x -> isequal(x, s), unit_scope)
            lab, nnd = cluster_units([last(items[u]) for u in us]; method = cfg.cluster_method,
                                     threshold = cfg.cluster_threshold, level = cfg.cluster_level,
                                     linkage = cfg.linkage, align_unequal = cfg.align_unequal,
                                     threaded = cfg.threaded, max_hierarchical = cfg.max_hierarchical)
            for (k, u) in enumerate(us)
                key = (s, lab[k])
                uclust[u] = get!(() -> newgroup!(sp === nothing ? lab[k] : key, s), gid, key)
                unn[u] = nnd[k]
            end
        end
        for r in 1:n
            u = unit[r]
            u === missing && continue
            row_cluster[r] = uclust[u]
            nn[r] = isfinite(unn[u]) ? unn[u] : missing
        end
    end
    # sizes in unique sequences and rows
    usets = Dict{Int,Set{Int}}(); rcount = Dict{Int,Int}()
    for r in 1:n
        c = row_cluster[r]
        c === missing && continue
        rcount[c] = get(rcount, c, 0) + 1
        unit[r] === missing || push!(get!(Set{Int}, usets, c), unit[r])
    end
    table = DataFrame(row = 1:n, id = ctx.ids,
                      species = sp === nothing ? fill(missing, n) : sp,
                      unit = unit, cluster = row_cluster,
                      cluster_label = Any[c === missing ? missing : labels[c] for c in row_cluster],
                      n_unique_in_cluster = Union{Missing,Int}[c === missing ? missing :
                                                               length(get(usets, c, Set{Int}())) for c in row_cluster],
                      n_rows_in_cluster = Union{Missing,Int}[c === missing ? missing : rcount[c] for c in row_cluster],
                      nn_distance = nn, cluster_source = fill(source, n),
                      cluster_basis = fill(basis, n),
                      threshold = fill(source in (:connected_components, :hierarchical) ?
                                       cfg.cluster_threshold : missing, n),
                      reason = row_reason)
    Groups(row_cluster, labels, gspecies, row_reason, source, table)
end

function _units!(unit, keys)
    idx = Dict{Any,Int}(); items = Any[]
    for r in eachindex(keys)
        k = keys[r]
        k === nothing && continue
        unit[r] = get!(idx, k) do
            push!(items, k); length(items)
        end
    end
    items
end

# --- Comparative analyses ----------------------------------------------------

"""Within-group results for one (gene, cluster) group."""
struct GroupResult
    gene::Symbol
    cluster::Int
    members::Vector{Int}            # all gene uniques of the group (frame resolved)
    comparable::Vector{Int}         # equal-length members used for comparisons
    status::Symbol                  # :ok or :insufficient
    reason::String
    metrics::Dict{Symbol,Vector{Union{Missing,Float64}}}   # aligned with `comparable`
    resolution::Dict{Symbol,Float64}                       # metric resolution floors
end

const COMPARATIVE_REASON_NONE = ""

"""
Outcome of the within-group comparison of one unique sequence: its amino-acid,
conservation and codon scores, or (`aa === nothing`) the reason why it was
not compared.
"""
struct GroupScores
    aa::Union{Nothing,NamedTuple}
    cons::Union{Nothing,NamedTuple}
    codon::Union{Nothing,NamedTuple}
    reason::String
end

GroupScores(reason::AbstractString) = GroupScores(nothing, nothing, nothing, String(reason))
scored(s::GroupScores) = s.aa !== nothing

const NOT_ANALYSED = GroupScores("not analysed")

function _mode_smallest(v::Vector{Int})
    c = Dict{Int,Int}()
    for x in v
        c[x] = get(c, x, 0) + 1
    end
    best = maximum(values(c))
    minimum(k for (k, n) in c if n == best)
end

function comparative_analyses!(ctx::AnalysisContext, L::Ledger, allgroups::Dict{Symbol,Groups})
    cfg = ctx.cfg
    aa_parts = DataFrame[]; cons_parts = DataFrame[]; cod_parts = DataFrame[]; grp_rows = NamedTuple[]
    results = GroupResult[]
    for g in ctx.genes
        groups = allgroups[genename(g)]
        res = Dict{Tuple{Int,Int},GroupScores}()       # (unique, cluster) => outcome
        bycluster = Dict{Int,Vector{Int}}()
        rows_of = Dict{Tuple{Int,Int},Vector{Int}}()
        for r in eachindex(g.row_unique)
            u = g.row_unique[r]; c = groups.row_cluster[r]
            (u == 0 || c === missing || g.records[u] === missing) && continue
            push!(get!(Vector{Int}, bycluster, c), u)
            push!(get!(Vector{Int}, rows_of, (u, c)), r)
        end
        for c in sort!(collect(keys(bycluster)))
            M = sort!(unique(bycluster[c]))
            gr, grow = _analyse_group(ctx, L, g, c, M, rows_of, res)
            push!(results, gr); push!(grp_rows, grow)
        end
        push!(aa_parts, _scores_rows(ctx, groups, g, res, :aa))
        push!(cons_parts, _scores_rows(ctx, groups, g, res, :cons))
        push!(cod_parts, _scores_rows(ctx, groups, g, res, :codon))
    end
    (aa_table = reduce(vcat, aa_parts; cols = :union),
     conservation_table = reduce(vcat, cons_parts; cols = :union),
     codon_table = reduce(vcat, cod_parts; cols = :union),
     group_table = isempty(grp_rows) ? DataFrame() :
         insertcols!(DataFrame(grp_rows), 2,
                     :species => [allgroups[x.gene].group_species[x.cluster] for x in grp_rows],
                     :cluster_label => [allgroups[x.gene].labels[x.cluster] for x in grp_rows]),
     results = results)
end

function _analyse_group(ctx, L, g::GeneData, c::Int, M::Vector{Int}, rows_of, res)
    cfg = ctx.cfg; gn = genename(g); mg = cfg.min_group_size
    recs = [g.records[u] for u in M]
    rowsfor(u) = get(rows_of, (u, c), Int[])
    nlen = [r.length_nt for r in recs]
    empty_metrics = Dict{Symbol,Vector{Union{Missing,Float64}}}()
    grow(status, reason; ncomp = 0, extra = NamedTuple()) = merge(
        (gene = gn, cluster = c, n_members = length(M), n_comparable = ncomp,
         status = status, reason = reason),
        merge((modal_length_nt = missing, n_codons = missing, n_effective = missing,
               var_frac_pos1 = missing, var_frac_pos2 = missing, var_frac_pos3 = missing,
               aa_variable_frac = missing, mean_conservation = missing), extra))
    if length(M) < mg
        reason = "group has $(length(M)) unique $(gn) sequence(s); at least $mg needed"
        for u in M
            res[(u, c)] = GroupScores(reason)
            for r in rowsfor(u)
                addflag!(L, r, gn, :insufficient_comparators; metric = "n_unique_in_group",
                         value = length(M), threshold = mg, detail = reason)
            end
        end
        return GroupResult(gn, c, M, Int[], :insufficient, reason, empty_metrics, Dict{Symbol,Float64}()),
               grow(:insufficient, reason)
    end
    # Length differences are assessed by the length test (length.jl).
    modal_nt = _mode_smallest(nlen)
    modal_cols = _mode_smallest([r.n_columns for r in recs])
    C = [u for (u, r) in zip(M, recs) if r.n_columns == modal_cols]
    for u in M
        g.records[u].n_columns == modal_cols && continue
        reason = "length ($(g.records[u].n_columns) columns) differs from the group modal length " *
                 "($modal_cols); positions are not compared without an alignment"
        res[(u, c)] = GroupScores(reason)
        for r in rowsfor(u)
            addflag!(L, r, gn, :not_comparable; detail = reason)
        end
    end
    if length(C) < mg
        reason = "only $(length(C)) equal-length $(gn) sequences in the group; at least $mg needed"
        for u in C
            res[(u, c)] = GroupScores(reason)
            for r in rowsfor(u)
                addflag!(L, r, gn, :insufficient_comparators; metric = "n_comparable_in_group",
                         value = length(C), threshold = mg, detail = reason)
            end
        end
        return GroupResult(gn, c, M, Int[], :insufficient, reason, empty_metrics, Dict{Symbol,Float64}()),
               grow(:insufficient, reason; ncomp = length(C))
    end
    # --- comparisons within the comparable set ---
    code = genetic_code(cfg.genetic_code); tabs = codon_tables(code)
    prots = [Vector{UInt8}(codeunits(g.records[u].protein)) for u in C]
    w = if cfg.conservation_weights === :henikoff
        henikoff_weights(prots)
    elseif cfg.conservation_weights === :abundance
        ctx.abundance === nothing && throw(ArgumentError("conservation_weights = :abundance needs abundance_col"))
        [sum(coalesce(ctx.abundance[r], 0.0) for r in rowsfor(u); init = 0.0) for u in C]
    else
        ones(length(C))
    end
    core = conservation_core(prots, w)
    frame = g.frame
    ncod = g.records[C[1]].n_codons
    cidx = [begin
                m = g.masks[u]; lead = frame - 1
                [codon_index(m[lead + 3k - 2], m[lead + 3k - 1], m[lead + 3k]) for k in 1:ncod]
            end for u in C]
    CW = zeros(64, ncod); CR = zeros(Int, ncod)
    for (ci, wi) in zip(cidx, w), k in 1:ncod
        x = ci[k]
        x == 0 && continue
        CW[x, k] += wi; CR[k] += 1
    end
    min_inf = mg - 1
    metrics = Dict{Symbol,Vector{Union{Missing,Float64}}}(
        m => Vector{Union{Missing,Float64}}(missing, length(C))
        for m in (:mean_kd, :frac_windows_above, :grantham_sum, :cw_grantham_sum,
                  :n_radical_invariant, :diff_pos12, :gc3, :nonsyn_diffs))
    for (i, u) in enumerate(C)
        aa, cs = loo_protein_scores(core, prots[i], w[i], min_inf, cfg.rare_threshold, cfg.radical)
        cc = CodonCompare()
        for k in 1:ncod
            own = cidx[i][k]
            own == 0 && continue
            CR[k] - 1 < min_inf && continue
            # Tolerances relative to the leave-one-out weight (as `_loo_position`)
            tot = sum(view(CW, :, k)) - w[i]
            tot > 0 || continue
            zero_tol = 1e-12 * tot; tie_tol = 1e-9 * tot
            best = 0; bestv = -1.0; tie = false
            for x in 1:64
                v = CW[x, k] - (x == own ? w[i] : 0.0)
                v < zero_tol && continue
                if v > bestv + tie_tol
                    best = x; bestv = v; tie = false
                elseif abs(v - bestv) <= tie_tol
                    tie = true
                end
            end
            (best == 0 || tie) && continue
            add_codon_pair!(cc, tabs, best, own)
        end
        cod = codon_summary(cc; dnds = cfg.estimate_dnds)
        res[(u, c)] = GroupScores(aa, merge(cs, (n_comparators = length(C) - 1,
                                                 n_effective_group = n_effective(w))),
                                  cod, "")
        h = g.hydro[u]
        metrics[:mean_kd][i] = h.mean_kd
        metrics[:frac_windows_above][i] = h.frac_windows_above
        metrics[:grantham_sum][i] = aa.grantham_sum
        metrics[:cw_grantham_sum][i] = cs.cw_grantham_sum
        metrics[:n_radical_invariant][i] = cs.n_radical_invariant
        metrics[:diff_pos12][i] = cc.n_codons == 0 ? missing : cod.diff_pos12
        metrics[:gc3][i] = g.records[u].gc3
        metrics[:nonsyn_diffs][i] = cc.n_codons == 0 ? missing : cod.nonsyn_diffs
    end
    # resolution floors: the smallest meaningful change of each metric
    nres = median([g.hydro[u].n_residues_scored for u in C])
    nwin = median([g.hydro[u].n_windows for u in C])
    resolution = Dict{Symbol,Float64}(
        :mean_kd => nres > 0 ? 4.5 / nres : 0.0,
        :frac_windows_above => nwin > 0 ? 1 / nwin : 0.0,
        :grantham_sum => 100.0, :cw_grantham_sum => 100.0, :n_radical_invariant => 1.0,
        :diff_pos12 => 1.0, :nonsyn_diffs => 1.0,
        :gc3 => ncod > 0 ? 1 / ncod : 0.0)
    # group-level variability
    var = zeros(Int, 3)
    for k in 1:ncod, j in 1:3
        seen = 0x00
        for u in C
            b = g.masks[u][frame - 1 + 3(k - 1) + j]
            is_unambiguous(b) && (seen |= b)
        end
        count_ones(seen) > 1 && (var[j] += 1)
    end
    aavar = count(p -> count(>(0), @view core.R[p, :]) > 1, 1:core.L)
    allc = [(_loo_position(core, p, 0, 0.0, mg)) for p in 1:core.L]
    cvals = [x[4] for x in allc if x[6]]
    extra = (modal_length_nt = modal_nt, n_codons = ncod, n_effective = n_effective(w),
             var_frac_pos1 = ncod > 0 ? var[1] / ncod : missing,
             var_frac_pos2 = ncod > 0 ? var[2] / ncod : missing,
             var_frac_pos3 = ncod > 0 ? var[3] / ncod : missing,
             aa_variable_frac = core.L > 0 ? aavar / core.L : missing,
             mean_conservation = isempty(cvals) ? missing : mean(cvals))
    GroupResult(gn, c, M, C, :ok, "", metrics, resolution), grow(:ok, ""; ncomp = length(C), extra)
end

const _AA_COLS = (:n_compared, :n_substitutions, :prop_different, :grantham_sum, :grantham_mean,
                  :grantham_max, :blosum62_mean, :blosum62_min, :delta_kd_sum, :abs_delta_kd_sum,
                  :n_radical, :n_stop_vs_consensus, :n_gap_vs_consensus, :n_positions_without_consensus)
const _CONS_COLS = (:n_scored_positions, :mean_conservation, :n_unseen_residues, :n_rare_residues,
                    :cw_grantham_sum, :cw_grantham_max, :cw_abs_delta_kd_sum, :n_radical_invariant,
                    :n_comparators, :n_effective_group)
const _CODON_COLS = (:n_codons_compared, :diff_pos1, :diff_pos2, :diff_pos3, :diff_pos12,
                     :n_diff_codons, :n_multi_diff_codons, :syn_diffs, :nonsyn_diffs, :prop_nonsyn,
                     :n_unclassified_codons, :syn_sites, :nonsyn_sites, :pN, :pS, :dN, :dS, :dN_dS,
                     :dnds_reason)

function _scores_rows(ctx, groups, g::GeneData, res, which::Symbol)
    n = length(g.row_unique)
    cols = which == :aa ? _AA_COLS : which == :cons ? _CONS_COLS : _CODON_COLS
    data = Dict(c => Vector{Any}(missing, n) for c in cols)
    reason = Vector{String}(undef, n)
    for r in 1:n
        u = g.row_unique[r]; c = groups.row_cluster[r]
        if u == 0
            reason[r] = string(g.row_status[r]); continue
        elseif g.records[u] === missing
            reason[r] = "reading frame unresolved"; continue
        elseif c === missing
            reason[r] = groups.row_reason[r] == "" ? "no comparison group" : groups.row_reason[r]; continue
        end
        v = get(res, (u, c), NOT_ANALYSED)
        if !scored(v)
            reason[r] = v.reason; continue
        end
        nt = getfield(v, which)
        for col in cols
            data[col][r] = getfield(nt, col)
        end
        reason[r] = ""
    end
    df = DataFrame(row = 1:n, id = ctx.ids, gene = fill(genename(g), n),
                   cluster = groups.row_cluster)
    for col in cols
        df[!, col] = _narrow(data[col])
    end
    df.reason = reason
    df
end

function _narrow(v::Vector{Any})
    T = mapreduce(typeof, promote_type, v; init = Missing)
    T = T === Missing ? Missing : T
    convert(Vector{T}, v)
end

"""Dataset codon usage per gene over unique gene sequences (frame resolved)."""
function codon_usage_table(gds::Vector{GeneData}, cfg::AssessmentConfig)
    code = genetic_code(cfg.genetic_code)
    parts = DataFrame[]
    for g in gds
        g.frame === nothing && continue
        counts = zeros(Int, 64)
        for u in 1:nunique(g)
            _count_codons!(counts, g.masks[u], g.frame)
        end
        df = _usage_table(counts, code)
        insertcols!(df, 1, :gene => fill(genename(g), 64))
        push!(parts, df)
    end
    isempty(parts) ? DataFrame() : reduce(vcat, parts)
end
