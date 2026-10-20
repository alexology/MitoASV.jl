# Analysis configuration. Every setting is recorded in `results.parameters`.

"""
    AssessmentConfig

All settings of an `assess_asvs` run. Built from the keyword arguments of
[`assess_asvs`](@ref); see its docstring for the meaning of each field.
"""
Base.@kwdef struct AssessmentConfig
    genetic_code::Int = 5
    terminal_stop_allowed::Bool = true
    frame_consensus_min::Float64 = 0.75
    frame_min_support::Int = 3
    # hydrophobicity
    hydrophobicity_window::Int = 19
    window_threshold::Float64 = 1.6
    hydrophobic_cutoff::Float64 = 0.0
    # grouping
    clustering::Bool = true
    cluster_method::Symbol = :connected_components
    cluster_threshold::Float64 = 0.03
    cluster_level::Symbol = :nucleotide
    linkage::Symbol = :average
    align_unequal::Bool = true
    cluster_on::Union{Nothing,Symbol} = nothing
    max_hierarchical::Int = 3000
    # length
    length_tolerance::Int = 0
    length_min_comparators::Int = 3
    # comparative analyses
    comparative::Bool = true
    min_group_size::Int = 5
    radical::RadicalCriteria = RadicalCriteria()
    conservation_weights::Union{Nothing,Symbol} = nothing
    rare_threshold::Float64 = 0.05
    estimate_dnds::Bool = false
    # anomalies
    anomaly_threshold::Dict{Symbol,Float64} = Dict{Symbol,Float64}()
    default_anomaly_threshold::Float64 = 3.5
    loo_dispersion_max::Int = 200
    # concordance
    concordance::Bool = true
    concordance_genes::Union{Nothing,Tuple{Symbol,Symbol}} = nothing
    concordance_scope::Symbol = :species
    concordance_k::Int = 5
    chimera_min_improvement::Float64 = 0.005
    pair_discordance_threshold::Float64 = 0.5
    max_concordance_asvs::Int = 3000
    max_discordant_pairs::Int = 1000
    compare_gene_trees::Bool = false
    # execution
    threaded::Bool = true
end

anomaly_threshold(cfg::AssessmentConfig, gene::Symbol) =
    get(cfg.anomaly_threshold, gene, cfg.default_anomaly_threshold)

function _validate(cfg::AssessmentConfig)
    genetic_code(cfg.genetic_code)
    0 < cfg.frame_consensus_min <= 1 || throw(ArgumentError("frame_consensus_min must be in (0, 1]"))
    cfg.hydrophobicity_window >= 1 || throw(ArgumentError("hydrophobicity_window must be ≥ 1"))
    cfg.cluster_method in (:connected_components, :hierarchical, :none) ||
        throw(ArgumentError("cluster_method must be :connected_components, :hierarchical or :none"))
    cfg.cluster_level in (:nucleotide, :protein) ||
        throw(ArgumentError("cluster_level must be :nucleotide or :protein"))
    cfg.linkage in (:average, :complete, :single) ||
        throw(ArgumentError("linkage must be :average, :complete or :single"))
    0 <= cfg.cluster_threshold < 1 || throw(ArgumentError("cluster_threshold must be in [0, 1)"))
    cfg.min_group_size >= 3 || throw(ArgumentError("min_group_size must be ≥ 3"))
    cfg.conservation_weights in (nothing, :henikoff, :abundance) ||
        throw(ArgumentError("conservation_weights must be nothing, :henikoff or :abundance"))
    cfg.concordance_scope in (:species, :all) ||
        throw(ArgumentError("concordance_scope must be :species or :all"))
    cfg.length_tolerance >= 0 || throw(ArgumentError("length_tolerance must be ≥ 0"))
    cfg.length_min_comparators >= 2 || throw(ArgumentError("length_min_comparators must be ≥ 2"))
    cfg.concordance_k >= 1 || throw(ArgumentError("concordance_k must be ≥ 1"))
    cfg.default_anomaly_threshold > 0 || throw(ArgumentError("anomaly threshold must be > 0"))
    cfg
end

function parameters_table(cfg::AssessmentConfig, extra::Vector{Pair{String,String}})
    names = String[]; values = String[]
    for f in fieldnames(AssessmentConfig)
        push!(names, String(f)); push!(values, _fmt(getfield(cfg, f)))
    end
    for (k, v) in extra
        push!(names, k); push!(values, v)
    end
    DataFrame(parameter = names, value = values)
end

_fmt(x) = string(x)
_fmt(::Nothing) = "nothing"
_fmt(d::Dict) = isempty(d) ? "default for all genes" :
    join(["$(k) => $(v)" for (k, v) in sort!(collect(d); by = first)], ", ")
