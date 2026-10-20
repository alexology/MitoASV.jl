# Robust within-group anomaly detection.

"""
    robust_anomaly_scores(values; resolution = 0.0, loo = true, min_comparators = 4)
        -> NamedTuple of vectors

Robust anomaly scores `(x - median) / dispersion` for each value, computed
against the **other** values (leave-one-out) when `loo = true`.

Dispersion, in order of preference:
1. `1.4826 × MAD` (consistent with the SD for normal data);
2. if the MAD is zero, `1.2533 × mean absolute deviation`
   (Iglewicz & Hoaglin 1993);
3. never below `resolution`, the smallest meaningful change of the metric
   (for example one substitution). This prevents extreme scores caused
   only by zero dispersion.
If the dispersion is still zero, the score is `0` for a value equal to the
median and `missing` otherwise (`reason` explains).

`missing` values, and values with fewer than `min_comparators` non-missing
comparators, get `missing` scores with a reason.

Returns `value`, `median`, `dispersion`, `dispersion_method`, `score`,
`reliable`, `n_comparators`, `reason`.
"""
function robust_anomaly_scores(values::AbstractVector; resolution::Real = 0.0, loo::Bool = true,
                               min_comparators::Integer = 4)
    n = length(values)
    med = Vector{Union{Missing,Float64}}(missing, n)
    disp = Vector{Union{Missing,Float64}}(missing, n)
    meth = fill("", n)
    score = Vector{Union{Missing,Float64}}(missing, n)
    reliable = falses(n)
    ncomp = zeros(Int, n)
    reason = fill("", n)
    allv = Float64[x for x in values if x !== missing && !isnan(x)]
    function centre(comps)
        m = median(comps)
        dev = abs.(comps .- m)
        d = 1.4826 * median(dev); method = "mad"
        if d == 0
            d = 1.2533 * mean(dev); method = "mean_absolute_deviation"
        end
        if d < resolution
            d = Float64(resolution); method = "resolution_floor"
        end
        (m, d, method)
    end
    global_centre = (!loo && length(allv) >= 1) ? centre(allv) : nothing
    for i in 1:n
        x = values[i]
        comps = loo ? Float64[values[j] for j in 1:n if j != i && values[j] !== missing && !isnan(values[j])] : allv
        ncomp[i] = length(comps)
        if x === missing || isnan(x)
            reason[i] = "value not available"; continue
        end
        if length(comps) < min_comparators
            reason[i] = "fewer than $min_comparators comparators with a value"; continue
        end
        m, d, method = loo ? centre(comps) : global_centre
        med[i] = m; disp[i] = d; meth[i] = method
        if d == 0
            if x == m
                score[i] = 0.0; reliable[i] = true; meth[i] = "zero_dispersion"
            else
                reason[i] = "zero dispersion among comparators; score not defined"
            end
            continue
        end
        score[i] = (x - m) / d
        reliable[i] = true
    end
    (value = Union{Missing,Float64}[x === missing ? missing : Float64(x) for x in values],
     median = med, dispersion = disp, dispersion_method = meth, score = score,
     reliable = collect(reliable), n_comparators = ncomp, reason = reason)
end

"""
Metrics scored within groups: `(metric, direction, flag, test)`.
`:upper` metrics only flag large values (more change than comparators).
"""
const ANOMALY_METRICS = [
    (:mean_kd, :two_sided, :hydrophobicity_outlier, :hydrophobicity),
    (:frac_windows_above, :two_sided, :hydrophobicity_outlier, :hydrophobicity),
    (:grantham_sum, :upper, :grantham_outlier, :grantham),
    (:cw_grantham_sum, :upper, :grantham_outlier, :grantham),
    (:diff_pos12, :upper, :codon_position_outlier, :codon_position),
    (:gc3, :two_sided, :codon_position_outlier, :codon_position),
    (:nonsyn_diffs, :upper, :nonsynonymous_outlier, :nonsynonymous),
]
const COMPARATIVE_TESTS = (:hydrophobicity, :grantham, :codon_position, :nonsynonymous)

_is_outlier(score, thr, dir) = score === missing ? missing :
    dir == :upper ? score >= thr : dir == :lower ? score <= -thr : abs(score) >= thr

function empty_anomaly_table()
    DataFrame(row = Int[], id = Any[], gene = Symbol[], cluster = Int[], metric = Symbol[],
              test = Symbol[], method = String[], direction = Symbol[],
              value = Union{Missing,Float64}[], group_median = Union{Missing,Float64}[],
              group_dispersion = Union{Missing,Float64}[], dispersion_method = String[],
              robust_score = Union{Missing,Float64}[], threshold = Float64[],
              is_outlier = Union{Missing,Bool}[], reliable = Bool[], n_comparators = Int[],
              reason = String[])
end

function anomaly_analysis!(ctx::AnalysisContext, L::Ledger, allgroups::Dict{Symbol,Groups}, comp)
    cfg = ctx.cfg
    buf = Tuple[]           # rows, converted to typed columns at the end
    gene_index = Dict(genename(g) => g for g in ctx.genes)
    evaluated = Set{Tuple{Int,Symbol,Symbol}}()      # (row, gene, test)
    for gr in comp.results
        gr.status == :ok || continue
        g = gene_index[gr.gene]
        thr = anomaly_threshold(cfg, gr.gene)
        rows_of = Dict{Int,Vector{Int}}()
        for r in eachindex(g.row_unique)
            u = g.row_unique[r]
            (u != 0 && allgroups[gr.gene].row_cluster[r] === gr.cluster) && push!(get!(Vector{Int}, rows_of, u), r)
        end
        loo = length(gr.comparable) <= cfg.loo_dispersion_max
        for (metric, dir, flag, test) in ANOMALY_METRICS
            s = robust_anomaly_scores(gr.metrics[metric]; resolution = gr.resolution[metric], loo,
                                      min_comparators = cfg.min_group_size - 1)
            for (i, u) in enumerate(gr.comparable)
                outl = _is_outlier(s.score[i], thr, dir)
                for r in get(rows_of, u, Int[])
                    push!(buf, (r, ctx.ids[r], gr.gene, gr.cluster, metric, test,
                                loo ? "robust_score_loo" : "robust_score", dir, s.value[i],
                                s.median[i], s.dispersion[i], s.dispersion_method[i], s.score[i],
                                thr, outl, s.reliable[i], s.n_comparators[i], s.reason[i]))
                    s.reliable[i] && push!(evaluated, (r, gr.gene, test))
                    if outl === true
                        addflag!(L, r, gr.gene, flag; metric = string(metric), value = s.score[i],
                                 threshold = thr,
                                 detail = "robust score $(round(s.score[i]; digits = 2)) " *
                                          "(value $(_r(s.value[i])), group median $(_r(s.median[i])), " *
                                          "dispersion $(_r(s.dispersion[i])) [$(s.dispersion_method[i])], $dir)")
                    end
                end
            end
        end
        # Rule: radical substitution at a site invariant among the comparators.
        vals = gr.metrics[:n_radical_invariant]
        for (i, u) in enumerate(gr.comparable), r in get(rows_of, u, Int[])
            v = vals[i]
            outl = v === missing ? missing : v >= 1
            push!(buf, (r, ctx.ids[r], gr.gene, gr.cluster, :n_radical_invariant, :grantham, "rule",
                        :upper, v, missing, missing, "", missing, 1.0, outl, v !== missing,
                        length(gr.comparable) - 1,
                        v === missing ? "value not available" : ""))
            if outl === true
                addflag!(L, r, gr.gene, :grantham_outlier; metric = "n_radical_invariant", value = v,
                         threshold = 1,
                         detail = "$(Int(v)) radical substitution(s) ($(cfg.radical)) at site(s) " *
                                  "invariant among the $(length(gr.comparable) - 1) comparators")
            end
        end
    end
    # Test bookkeeping: a comparative test is evaluated only when at least one
    # of its metrics was reliably scored for that row and gene.
    aa = comp.aa_table
    for k in 1:nrow(aa)
        r = aa.row[k]; gn = aa.gene[k]
        for test in COMPARATIVE_TESTS
            if (r, gn, test) in evaluated
                addtest!(L, r, gn, test, :evaluated)
            else
                reason = aa.reason[k] == "" ? "no reliable within-group score" : aa.reason[k]
                addtest!(L, r, gn, test, :not_evaluated; reason)
            end
        end
    end
    out = empty_anomaly_table()
    for (k, nm) in enumerate(propertynames(out))
        append!(out[!, nm], (row[k] for row in buf))
    end
    sort!(out, [:row, :gene, :metric])
end

_r(x) = x === missing ? "missing" : string(round(x; sigdigits = 4))
