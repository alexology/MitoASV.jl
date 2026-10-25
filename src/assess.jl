# The DataFrame interface: assess_asvs and assess_asvs_df.

"""
    AssessmentResult

Output of [`assess_asvs`](@ref). All tables are `DataFrame`s; per-ASV tables
have one row per input row (and gene), keyed by `row` (the input row index)
and `id`.

| Field | Content |
|---|---|
| `summary` | one row per input row: descriptive classification and flag lists |
| `validation` | reading frame, stops, ambiguous codons, length, GC by codon position (row × gene) |
| `composition` | amino acid composition (row × gene) |
| `hydrophobicity` | Kyte-Doolittle metrics (row × gene) |
| `lengths` | gene length against the length inferred from the other sequences of the species (row × gene) |
| `amino_acid_changes` | substitutions against the leave-one-out group consensus (row × gene) |
| `conservation` | leave-one-out conservation scores (row × gene) |
| `codon_statistics` | codon-aware comparison with the group consensus (row × gene) |
| `codon_usage` | dataset codon usage and RSCU (gene × codon) |
| `clusters` | comparison group of each row and gene |
| `group_statistics` | variability per group and gene |
| `gene_concordance` | per-row two-gene concordance (multi-gene runs) |
| `discordant_pairs` | ASV pairs with discordant gene distances |
| `anomalies` | robust within-group scores (row × gene × metric) |
| `flags` | all flags in long format |
| `tests` | which tests were evaluated, not evaluated (with reason) or not applicable |
| `parameters` | every setting used |
| `input` | the input DataFrame (not modified) |
"""
struct AssessmentResult
    summary::DataFrame
    validation::DataFrame
    composition::DataFrame
    hydrophobicity::DataFrame
    lengths::DataFrame
    amino_acid_changes::DataFrame
    conservation::DataFrame
    codon_statistics::DataFrame
    codon_usage::DataFrame
    clusters::DataFrame
    group_statistics::DataFrame
    gene_concordance::DataFrame
    discordant_pairs::DataFrame
    anomalies::DataFrame
    flags::DataFrame
    tests::DataFrame
    parameters::DataFrame
    input::DataFrame
end

function Base.show(io::IO, ::MIME"text/plain", r::AssessmentResult)
    n = nrow(r.summary)
    println(io, "MitoASV assessment of $n input rows")
    for c in (CLASS_NONE, CLASS_ONE, CLASS_MULTI, CLASS_INSUFFICIENT)
        println(io, "  ", rpad(c * ":", 34), count(==(c), r.summary.classification))
    end
    print(io, "Tables: summary, validation, composition, hydrophobicity, lengths, amino_acid_changes, ",
          "conservation, codon_statistics, codon_usage, clusters, group_statistics, ",
          "gene_concordance, discordant_pairs, anomalies, flags, tests, parameters")
end

function _gene_specs(gene, frame, strand, genes)
    if genes === nothing
        return [GeneSpec(gene === nothing ? :gene : Symbol(gene); frame, strand)]
    end
    gene === nothing || gene == :gene ||
        @warn "both `gene` and `genes` given; `gene` is ignored"
    itr = genes isa AbstractDict ? collect(genes) : genes isa AbstractVector ? genes : collect(pairs(genes))
    # `frame` and `strand` given as arguments apply to genes without their own
    # (a `GeneSpec` always has a strand, so only columns and NamedTuples inherit it)
    specs = [GeneSpec(first(p), _inherit_strand(last(p), strand)) for p in itr]
    frame === nothing || (specs = [s.frame === nothing ? GeneSpec(s.name, s.start, s.stop, frame, s.strand, s.column) : s
                                   for s in specs])
    isempty(specs) && throw(ArgumentError("`genes` is empty"))
    length(unique(s.name for s in specs)) == length(specs) ||
        throw(ArgumentError("gene names must be unique"))
    # Order: AbstractDict keys have no reliable order, so sort by position.
    genes isa AbstractDict && sort!(specs; by = s -> (coalesce(string(s.column), ""), s.start, string(s.name)))
    specs
end

_inherit_strand(column::Union{Symbol,AbstractString}, strand) = (column = Symbol(column), strand)
_inherit_strand(nt::NamedTuple, strand) = haskey(nt, :strand) ? nt : merge(nt, (; strand))
_inherit_strand(g, _) = g

_coalesce_nothing(x) = x === nothing ? missing : x

"""
    assess_asvs(df; sequence_col = :sequence, id_col = nothing, gene = :gene,
                genes = nothing, frame = nothing, strand = :+, genetic_code = 5,
                abundance_col = nothing, sample_col = nothing, replicate_col = nothing,
                cluster_col = nothing, species_col = nothing, kwargs...) -> AssessmentResult

Reference-free biological plausibility assessment of protein-coding ASVs.
The input `DataFrame` is not modified and no row is removed; every output
table maps back to input rows through `row`.

# Input
- **One column per gene** (the usual input): `genes = Dict(:ND2 => :nd2_seq,
  :COI => :coi_seq)` maps each gene to the column holding its sequences
  (strings, `BioSequences` sequences or `missing`). Columns can also be
  given as `(column = :nd2_seq, frame = 1, strand = :+)`. `frame` and
  `strand` passed as arguments apply to genes without their own. Without a
  full-sequence column, each gene is grouped on its own sequences, so a row
  missing one gene is still analysed for the others.
- `sequence_col`: column with the whole amplicon, for single-gene analyses
  or when genes are given by coordinates. It is ignored if absent and
  every gene has a column.
- `id_col`: optional ASV identifiers (default: the row number).
- `gene`, `frame`, `strand`: single-gene analysis of the whole sequence.
- `genes` with coordinates on `sequence_col`, e.g.
  `Dict(:ND2 => (start = 1, stop = 900, frame = 1), :COI => (start = 1200, stop = 2700, frame = 1))`.
  Coordinates are 1-based and inclusive; see [`GeneSpec`](@ref).
- Length: each gene's length is compared with the expected length, inferred
  as the modal length of the other unique sequences of that gene in the
  same species (or comparison group): `length_tolerance = 0` nt,
  `length_min_comparators = 3`. See `results.lengths`.
- `genetic_code`: NCBI table (default 5, invertebrate mitochondrial).
- `abundance_col`: read counts; used for chimera parent abundance and, if
  requested, abundance weights. With `sample_col`, the same sequence may
  appear in several rows (long format): sequences are analysed once.
- `sample_col`, `replicate_col`: reported as detection support (number of
  samples/replicates), never used as evidence of authenticity.
- `species_col`: species of each ASV (for example a `species` column from a
  taxonomic classifier). When given, **every comparison stays within a
  species**: clustering is run separately in each species, so groups never
  mix species; with `clustering = false` each species is one group; with
  `cluster_col` the group is the (species, cluster) pair; the ND2–COI
  concordance compares ASVs of the same species (`concordance_scope =
  :species`; `:all` compares across species, needed to detect
  interspecific chimeras). A sequence assigned to two species is analysed
  in each separately. Rows with a missing species get no comparison group
  (`:missing_species`) and their comparative tests are not evaluated. The
  reading-frame consensus remains dataset-wide (it depends on the primers).
- `cluster_col`: user-provided comparison groups (replaces clustering;
  combined with `species_col` when both are given).

# Settings (keyword arguments, all recorded in `results.parameters`)
- Frames: `terminal_stop_allowed = true`, `frame_consensus_min = 0.75`,
  `frame_min_support = 3`.
- Hydrophobicity: `hydrophobicity_window = 19`, `window_threshold = 1.6`,
  `hydrophobic_cutoff = 0.0`.
- Grouping: `clustering = true`, `cluster_method = :connected_components`
  (or `:hierarchical`, `:none`), `cluster_threshold = 0.03`,
  `cluster_level = :nucleotide` (or `:protein`), `linkage = :average`,
  `align_unequal = true`, `cluster_on = nothing` (a gene name to cluster all
  genes on that gene, or `:per_gene`; default: the full sequence if there is
  a `sequence_col`, otherwise each gene separately), `max_hierarchical = 3000`.
- Comparative analyses: `comparative = true`, `min_group_size = 5`,
  `radical = RadicalCriteria()`, `conservation_weights = nothing`
  (`:henikoff` or `:abundance`), `rare_threshold = 0.05`,
  `estimate_dnds = false`.
- Anomalies: `anomaly_threshold = 3.5` (a number, or a `Dict` gene =>
  threshold for gene-specific calibration), `loo_dispersion_max = 200`.
- Concordance: `concordance = true`, `concordance_genes = nothing` (first
  two genes), `concordance_k = 5`, `chimera_min_improvement = 0.005`,
  `concordance_scope = :species`, `pair_discordance_threshold = 0.5`, `max_concordance_asvs = 3000`,
  `max_discordant_pairs = 1000`, `compare_gene_trees = false`.
- `threaded = true`: use Julia threads (results do not depend on the
  number of threads).
"""
function assess_asvs(df::AbstractDataFrame; sequence_col = :sequence, id_col = nothing,
                     gene = :gene, genes = nothing, frame = nothing, strand = :+,
                     genetic_code = 5, abundance_col = nothing, sample_col = nothing,
                     replicate_col = nothing, cluster_col = nothing, species_col = nothing,
                     anomaly_threshold = 3.5, kwargs...)
    input = DataFrame(df; copycols = false)
    n = nrow(input)
    athr = anomaly_threshold isa AbstractDict ?
        Dict{Symbol,Float64}(Symbol(k) => Float64(v) for (k, v) in anomaly_threshold) :
        Dict{Symbol,Float64}()
    default_thr = anomaly_threshold isa Real ? Float64(anomaly_threshold) : 3.5
    cfg = _validate(AssessmentConfig(; genetic_code = Int(genetic_code isa GeneticCode ?
                                                         genetic_code.id : genetic_code),
                                     anomaly_threshold = athr,
                                     default_anomaly_threshold = default_thr, kwargs...))
    specs = _gene_specs(gene, frame, strand, genes)
    for c in (id_col, abundance_col, sample_col, replicate_col, cluster_col, species_col)
        c === nothing || hasproperty(input, c) ||
            throw(ArgumentError("column $(repr(c)) not found in the input DataFrame"))
    end
    ids = id_col === nothing ? collect(1:n) : collect(input[!, id_col])
    needs_main = any(s -> s.column === nothing, specs) ||
                 (cfg.cluster_on === nothing && sequence_col !== nothing && hasproperty(input, sequence_col))
    main = (sequence_col !== nothing && (needs_main || any(s -> s.column === nothing, specs))) ?
           _parse_column(input, sequence_col) : nothing
    abundance = abundance_col === nothing ? nothing :
        Union{Missing,Float64}[x === missing ? missing : Float64(x) for x in input[!, abundance_col]]

    colcache = Dict{Symbol,Any}()
    gds = [prepare_gene(s, input, main, cfg, colcache) for s in specs]
    L = Ledger()

    _coding_flags!(L, gds, main, colcache, cfg)
    ctx = AnalysisContext(input, ids, cfg, gds, main, abundance, colcache)

    groups = assign_groups!(ctx, L, cluster_col, species_col)
    lengths = length_analysis!(ctx, L, groups, species_col)
    comp = cfg.comparative ? comparative_analyses!(ctx, L, groups) : nothing
    anomalies = cfg.comparative ? anomaly_analysis!(ctx, L, groups, comp) : empty_anomaly_table()
    if !cfg.comparative
        for g in gds, r in 1:n, t in COMPARATIVE_TESTS
            addtest!(L, r, genename(g), t, :not_evaluated; reason = "disabled (comparative = false)")
        end
    end
    conc = concordance_analysis!(ctx, L, species_col)

    summary = composite_summary(L, n, ids)
    species_col === nothing || insertcols!(summary, 3, :species => collect(input[!, species_col]))
    _add_support!(summary, input, ctx, sample_col, replicate_col, abundance)
    extra = ["genes" => join([_spec_string(s) for s in specs], "; "),
             "frames_used" => join(["$(genename(g)) => $(g.frame === nothing ? "unresolved" : g.frame) ($(g.frame_source))"
                                    for g in gds], "; "),
             "frame_notes" => join(["$(genename(g)): $(g.frame_note)" for g in gds], " | "),
             "sequence_col" => string(sequence_col), "id_col" => string(id_col),
             "abundance_col" => string(abundance_col), "sample_col" => string(sample_col),
             "replicate_col" => string(replicate_col), "cluster_col" => string(cluster_col),
             "species_col" => string(species_col),
             "threads" => string(cfg.threaded ? Threads.nthreads() : 1),
             "package_version" => string(pkgversion(MitoASV))]
    AssessmentResult(summary,
                     reduce(vcat, [validation_table(g, ids, _invalid_counts(g, main, colcache)) for g in gds]),
                     reduce(vcat, [composition_rows(g, ids) for g in gds]),
                     reduce(vcat, [hydrophobicity_rows(g, ids) for g in gds]),
                     lengths,
                     comp === nothing ? DataFrame() : comp.aa_table,
                     comp === nothing ? DataFrame() : comp.conservation_table,
                     comp === nothing ? DataFrame() : comp.codon_table,
                     codon_usage_table(gds, cfg),
                     clusters_table(groups, gds),
                     comp === nothing ? DataFrame() : comp.group_table,
                     conc.table, conc.pairs, anomalies,
                     flags_table(L, ids), tests_table(L, ids),
                     parameters_table(cfg, extra), input)
end

_spec_string(s::GeneSpec) = "$(s.name): " *
    (s.column === nothing ? "" : "column=$(s.column), ") *
    "start=$(s.start), stop=$(s.stop === nothing ? "end" : s.stop), " *
    "frame=$(s.frame === nothing ? "unknown" : s.frame), strand=$(s.strand)"

function _invalid_counts(g::GeneData, main, colcache)
    src = g.spec.column === nothing ? main : colcache[g.spec.column]
    Union{Missing,Int}[p === missing ? missing : p.n_invalid for p in src]
end


function _coding_flags!(L::Ledger, gds, main, colcache, cfg)
    for g in gds
        gn = genename(g)
        inv = _invalid_counts(g, main, colcache)
        for r in eachindex(g.row_unique)
            st = g.row_status[r]
            inv[r] !== missing && inv[r] > 0 &&
                addflag!(L, r, gn, :invalid_symbols; metric = "n_invalid_symbols", value = inv[r],
                         detail = "non-IUPAC symbols read as N; the input is unchanged")
            if st in (:missing_sequence, :empty_sequence, :absent)
                addflag!(L, r, gn, st == :absent ? :gene_absent : st;
                         detail = st == :absent ? "the sequence ends before the gene start" : "")
                addtest!(L, r, gn, :coding, :not_evaluated; reason = string(st))
                continue
            end
            st == :truncated && addflag!(L, r, gn, :gene_truncated;
                                         detail = "the sequence ends before the gene stop; analysed as a partial gene")
            u = g.row_unique[r]
            rec = g.records[u]
            if rec === missing
                addflag!(L, r, gn, :frame_unresolved; detail = g.frame_note)
                addtest!(L, r, gn, :coding, :not_evaluated; reason = "reading frame unresolved")
                continue
            end
            addtest!(L, r, gn, :coding, :evaluated)
            if rec.n_internal_stops > 0
                addflag!(L, r, gn, :internal_stop; metric = "n_internal_stops",
                         value = rec.n_internal_stops, threshold = 0,
                         detail = "stop codon(s) at codon position(s) " *
                                  join(rec.internal_stop_positions, ",") * " in frame $(rec.frame)")
            end
            if rec.n_codons == 0
                addflag!(L, r, gn, :invalid_frame; metric = "n_codons", value = 0,
                         detail = "no complete codon in frame $(rec.frame)")
            elseif rec.n_partial_gap_codons > 0
                addflag!(L, r, gn, :invalid_frame; metric = "n_partial_gap_codons",
                         value = rec.n_partial_gap_codons, threshold = 0,
                         detail = "codon(s) with 1 or 2 gap columns: frameshift-like indel or misalignment")
            end
            rec.n_possible_stops > 0 &&
                addflag!(L, r, gn, :possible_stop; metric = "n_possible_stops",
                         value = rec.n_possible_stops,
                         detail = "ambiguous codon(s) that may encode a stop")
        end
    end
end

"""
    assess_asvs_df(df; kwargs...) -> DataFrame

Run [`assess_asvs`](@ref) and return one row per input row: the input
columns, the composite summary, one Bool column per anomaly flag
(`missing` when the corresponding tests could not be evaluated and no flag
was raised) and the main metrics of every gene, prefixed with the gene name
(for example `COI_mean_kd`). Column names already used by the input get a
`mito_` prefix, repeated until the generated name is unused.
"""
assess_asvs_df(df::AbstractDataFrame; kwargs...) = to_wide(assess_asvs(df; kwargs...))
