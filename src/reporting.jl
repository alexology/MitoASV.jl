# Flags, test bookkeeping and the descriptive composite summary.

"""
Type and evidence class of every flag. Only `:anomaly` flags count towards
the composite summary; `:data` and `:status` flags describe the input or
why a test could not be run.
"""
const FLAG_INFO = Dict{Symbol,Tuple{Symbol,Union{Missing,Symbol}}}(
    :internal_stop => (:anomaly, :coding_integrity),
    :invalid_frame => (:anomaly, :coding_integrity),
    :length_outlier => (:anomaly, :coding_integrity),
    :hydrophobicity_outlier => (:anomaly, :biochemical),
    :grantham_outlier => (:anomaly, :biochemical),
    :codon_position_outlier => (:anomaly, :evolutionary),
    :nonsynonymous_outlier => (:anomaly, :evolutionary),
    :gene_discordance => (:anomaly, :gene_concordance),
    :insufficient_comparators => (:status, missing),
    :frame_unresolved => (:status, missing),
    :not_comparable => (:status, missing),
    :missing_sequence => (:data, missing),
    :missing_species => (:data, missing),
    :empty_sequence => (:data, missing),
    :invalid_symbols => (:data, missing),
    :gene_absent => (:data, missing),
    :gene_truncated => (:data, missing),
    :possible_stop => (:data, missing),
)

const ANOMALY_FLAGS = [:internal_stop, :invalid_frame, :length_outlier, :hydrophobicity_outlier, :grantham_outlier,
                       :codon_position_outlier, :nonsynonymous_outlier, :gene_discordance]

struct FlagRecord
    row::Int
    gene::Union{Missing,Symbol}
    flag::Symbol
    metric::String
    value::Union{Missing,Float64}
    threshold::Union{Missing,Float64}
    detail::String
end

struct TestRecord
    row::Int
    gene::Union{Missing,Symbol}
    test::Symbol
    status::Symbol          # :evaluated, :not_evaluated, :not_applicable
    reason::String
end

mutable struct Ledger
    flags::Vector{FlagRecord}
    tests::Vector{TestRecord}
end
Ledger() = Ledger(FlagRecord[], TestRecord[])

addflag!(L::Ledger, row, gene, flag; metric = "", value = missing, threshold = missing, detail = "") =
    push!(L.flags, FlagRecord(row, gene, flag, metric,
                              value === missing ? missing : Float64(value),
                              threshold === missing ? missing : Float64(threshold), detail))

addtest!(L::Ledger, row, gene, test, status; reason = "") =
    push!(L.tests, TestRecord(row, gene, test, status, reason))

_flag_order(f::Symbol) = (k = findfirst(==(f), ANOMALY_FLAGS); k === nothing ? 100 : k)

function flags_table(L::Ledger, ids)
    fl = sort(L.flags; by = f -> (f.row, coalesce(string(f.gene), ""), _flag_order(f.flag), string(f.flag), f.metric))
    DataFrame(row = [f.row for f in fl], id = [ids[f.row] for f in fl],
              gene = Union{Missing,Symbol}[f.gene for f in fl], flag = [f.flag for f in fl],
              flag_type = [FLAG_INFO[f.flag][1] for f in fl],
              evidence_class = Union{Missing,Symbol}[FLAG_INFO[f.flag][2] for f in fl],
              metric = [f.metric for f in fl], value = Union{Missing,Float64}[f.value for f in fl],
              threshold = Union{Missing,Float64}[f.threshold for f in fl],
              detail = [f.detail for f in fl])
end

function tests_table(L::Ledger, ids)
    t = sort(L.tests; by = x -> (x.row, coalesce(string(x.gene), ""), string(x.test)))
    DataFrame(row = [x.row for x in t], id = [ids[x.row] for x in t],
              gene = Union{Missing,Symbol}[x.gene for x in t], test = [x.test for x in t],
              status = [x.status for x in t], reason = [x.reason for x in t])
end

const CLASS_NONE = "No detected anomaly"
const CLASS_ONE = "One or more anomalies"
const CLASS_MULTI = "Multiple independent anomalies"
const CLASS_INSUFFICIENT = "Insufficient evidence"

"""
Descriptive composite summary per row. Heterogeneous metrics are never
summed. The composite is the number of distinct evidence classes with at
least one anomaly flag (coding integrity, biochemical, evolutionary, gene
concordance) — not a probability. Biochemical or evolutionary flags on a
gene that also has a coding-integrity flag (internal stop, frameshift-like
length) follow from that defect; they are listed in `dependent_flags` and
not counted as independent evidence. A test that was not evaluated is never
counted as passed: a row without anomaly flags is "No detected anomaly"
only when every applicable test was evaluated, otherwise
"Insufficient evidence".
"""
function composite_summary(L::Ledger, n::Int, ids)
    flags_by = [FlagRecord[] for _ in 1:n]
    for f in L.flags
        push!(flags_by[f.row], f)
    end
    tests_by = [TestRecord[] for _ in 1:n]
    for t in L.tests
        push!(tests_by[t.row], t)
    end
    label(gene, x) = gene === missing ? string(x) : "$(gene):$(x)"
    classification = String[]; n_flags = Int[]; n_classes = Int[]
    flag_list = String[]; class_list = String[]; status_list = String[]
    n_eval = Int[]; n_noteval = Int[]; noteval_list = String[]; dep_list = String[]
    for r in 1:n
        an = [f for f in flags_by[r] if FLAG_INFO[f.flag][1] == :anomaly]
        anames = unique([label(f.gene, f.flag) for f in an])
        # Comparative flags on a gene that has a coding-integrity defect are
        # consequences of that defect: listed, but not independent evidence.
        coding_genes = Set(f.gene for f in an if FLAG_INFO[f.flag][2] == :coding_integrity)
        isdep(f) = FLAG_INFO[f.flag][2] != :coding_integrity && f.gene !== missing && f.gene in coding_genes
        classes = sort(unique([string(FLAG_INFO[f.flag][2]) for f in an if !isdep(f)]))
        push!(dep_list, join(unique([label(f.gene, f.flag) for f in an if isdep(f)]), ";"))
        other = unique([label(f.gene, f.flag) for f in flags_by[r] if FLAG_INFO[f.flag][1] != :anomaly])
        ev = count(t -> t.status == :evaluated, tests_by[r])
        ne = [t for t in tests_by[r] if t.status == :not_evaluated]
        push!(n_flags, length(anames)); push!(n_classes, length(classes))
        push!(flag_list, join(anames, ";")); push!(class_list, join(classes, ";"))
        push!(status_list, join(other, ";"))
        push!(n_eval, ev); push!(n_noteval, length(ne))
        push!(noteval_list, join(unique([label(t.gene, t.test) for t in ne]), ";"))
        c = length(classes) >= 2 ? CLASS_MULTI :
            length(anames) >= 1 ? CLASS_ONE :
            (isempty(ne) && ev > 0) ? CLASS_NONE : CLASS_INSUFFICIENT
        push!(classification, c)
    end
    DataFrame(row = 1:n, id = ids, classification = classification,
              n_anomaly_flags = n_flags, n_evidence_classes = n_classes,
              anomaly_flags = flag_list, evidence_classes = class_list, dependent_flags = dep_list,
              n_tests_evaluated = n_eval, n_tests_not_evaluated = n_noteval,
              tests_not_evaluated = noteval_list, other_flags = status_list)
end

"""
Detection support of each row's sequence: rows sharing the sequence, total
abundance, distinct samples and replicates. Descriptive only; never used
as evidence of authenticity.
"""
function _add_support!(summary::DataFrame, input::DataFrame, ctx, sample_col, replicate_col, abundance)
    n = nrow(input)
    key = Vector{Any}(nothing, n)
    for r in 1:n
        if ctx.main !== nothing
            p = ctx.main[r]
            (p === missing || isempty(p.mask)) || (key[r] = p.mask)
        else
            all(g -> g.row_unique[r] != 0, ctx.genes) &&
                (key[r] = Tuple(g.row_unique[r] for g in ctx.genes))
        end
    end
    groups = Dict{Any,Vector{Int}}()
    for r in 1:n
        key[r] === nothing || push!(get!(Vector{Int}, groups, key[r]), r)
    end
    detected(r) = abundance === nothing || (abundance[r] !== missing && abundance[r] > 0)
    nrows = Vector{Union{Missing,Int}}(missing, n)
    tot = Vector{Union{Missing,Float64}}(missing, n)
    nsamp = Vector{Union{Missing,Int}}(missing, n)
    nrep = Vector{Union{Missing,Int}}(missing, n)
    for rows in values(groups)
        t = abundance === nothing ? missing :
            (any(r -> abundance[r] === missing, rows) ? missing : sum(abundance[r] for r in rows))
        det = filter(detected, rows)
        ns = sample_col === nothing ? missing : length(unique(input[r, sample_col] for r in det))
        nr = replicate_col === nothing ? missing :
             length(unique((sample_col === nothing ? nothing : input[r, sample_col],
                            input[r, replicate_col]) for r in det))
        for r in rows
            nrows[r] = length(rows); tot[r] = t; nsamp[r] = ns; nrep[r] = nr
        end
    end
    summary.n_rows_same_sequence = nrows
    summary.total_abundance = tot
    summary.n_samples_detected = nsamp
    summary.n_replicates_detected = nrep
    summary
end

const FLAG_TESTS = Dict(:internal_stop => :coding, :invalid_frame => :coding,
                        :length_outlier => :length,
                        :hydrophobicity_outlier => :hydrophobicity, :grantham_outlier => :grantham,
                        :codon_position_outlier => :codon_position,
                        :nonsynonymous_outlier => :nonsynonymous,
                        :gene_discordance => :gene_concordance)

const WIDE_COLUMNS = [
    (:lengths, [:modal_length_nt, :length_difference]),
    (:validation, [:status, :frame, :frame_source, :length_nt, :n_internal_stops, :terminal_stop,
                   :trailing_bases, :n_ambiguous_codons, :gc, :gc3]),
    (:hydrophobicity, [:mean_kd, :sd_kd, :frac_hydrophobic, :frac_windows_above, :max_window_kd]),
    (:amino_acid_changes, [:n_substitutions, :grantham_sum, :grantham_max, :n_radical]),
    (:conservation, [:n_unseen_residues, :cw_grantham_sum, :n_radical_invariant]),
    (:codon_statistics, [:diff_pos12, :diff_pos3, :syn_diffs, :nonsyn_diffs, :prop_nonsyn]),
]

"""One row per input row: input columns, summary, flag columns and gene metrics."""
function to_wide(res)
    n = nrow(res.input)
    out = copy(res.input)
    used = Set(propertynames(out))
    put!(name, v) = begin
        nm = name in used ? Symbol("mito_", name) : name
        out[!, nm] = v; push!(used, nm)
    end
    for c in propertynames(res.summary)
        c in (:row, :id, :species) && continue
        put!(c, res.summary[!, c])
    end
    for gdf in groupby(res.clusters, :gene; sort = true)
        put!(Symbol(gdf.gene[1], "_cluster"), gdf.cluster)
    end
    # flag columns: true if raised; false only if the underlying test was
    # evaluated for every applicable gene; otherwise missing.
    raised = Set{Tuple{Int,Symbol}}((r, f) for (r, f) in zip(res.flags.row, res.flags.flag))
    tests = Dict{Tuple{Int,Symbol},Vector{Symbol}}()
    for (r, t, s) in zip(res.tests.row, res.tests.test, res.tests.status)
        push!(get!(Vector{Symbol}, tests, (r, t)), s)
    end
    for f in ANOMALY_FLAGS
        t = FLAG_TESTS[f]
        v = Vector{Union{Missing,Bool}}(missing, n)
        for r in 1:n
            if (r, f) in raised
                v[r] = true
            else
                st = get(tests, (r, t), Symbol[])
                if !isempty(st) && all(==(:evaluated), st)
                    v[r] = false
                elseif !isempty(st) && all(==(:not_applicable), st)
                    v[r] = missing
                end
            end
        end
        put!(Symbol("flag_", f), v)
    end
    for (tbl, cols) in WIDE_COLUMNS
        t = getfield(res, tbl)
        nrow(t) == 0 && continue
        for gdf in groupby(t, :gene; sort = true)
            gname = gdf.gene[1]
            idx = Dict(r => k for (k, r) in enumerate(gdf.row))
            for c in cols
                hasproperty(gdf, c) || continue
                col = gdf[!, c]
                put!(Symbol(gname, "_", c),
                     [haskey(idx, r) ? col[idx[r]] : missing for r in 1:n])
            end
        end
    end
    gc = res.gene_concordance
    if hasproperty(gc, :two_parent_improvement)
        for c in (:two_parent_improvement, :knn_jaccard, :row_spearman_nt, :nn_gene1_id, :nn_gene2_id)
            put!(Symbol("concordance_", c), gc[!, c])
        end
    end
    out
end
