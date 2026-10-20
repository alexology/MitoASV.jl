"""
    MitoASV

Reference-free biological plausibility assessment of mitochondrial
protein-coding amplicon sequence variants (ASVs).

Plausibility is not authenticity: the package reports interpretable scores,
diagnostic flags and supporting statistics. It never decides that an ASV is
genuine or artifactual, and without external references, read-level data or
experimental evidence it cannot separate genuine variants from sequencing
errors, NUMTs or chimeras.
"""
module MitoASV

using DataFrames
using Statistics

export assess_asvs, assess_asvs_df, AssessmentResult, AssessmentConfig, GeneSpec,
       # sequences and translation
       GeneticCode, genetic_code, parse_sequence, translate_sequence, translate_codon,
       coding_validation, frame_candidates, aa_composition,
       # amino acid properties
       KYTE_DOOLITTLE, grantham_distance, blosum62_score, RadicalCriteria, is_radical,
       hydrophobicity_profile, hydrophobicity_metrics, hydrophobicity_table,
       compare_proteins,
       # codons, grouping, conservation
       compare_codons, codon_usage, cluster_asvs, conservation_analysis, henikoff_weights,
       # anomalies, concordance
       robust_anomaly_scores, gene_concordance, spearman, modal_length,
       simulate_asv_dataset

include("geneticcode.jl")
include("sequences.jl")
include("validation.jl")
include("aaproperties.jl")
include("hydrophobicity.jl")
include("substitutions.jl")
include("codons.jl")
include("distances.jl")
include("clustering.jl")
include("conservation.jl")
include("config.jl")
include("reporting.jl")
include("pipeline_genes.jl")
include("pipeline_groups.jl")
include("length.jl")
include("anomaly.jl")
include("concordance.jl")
include("assess.jl")
include("simulate.jl")

end
