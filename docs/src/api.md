# API

```@meta
CurrentModule = MitoASV
```

```@docs
MitoASV
```

## Assessment

```@docs
AssessmentResult
assess_asvs
assess_asvs_df
```

## Settings

```@docs
AssessmentConfig
```

## Genetic codes

```@docs
GeneticCode
genetic_code
translate_codon
```

## Sequences and genes

```@docs
GeneSpec
parse_sequence
```

## Translation and coding validation

```@docs
aa_composition
coding_validation
frame_candidates
translate_sequence
```

## Amino-acid properties

```@docs
KYTE_DOOLITTLE
RadicalCriteria
blosum62_score
grantham_distance
is_radical
```

## Hydrophobicity

```@docs
hydrophobicity_metrics
hydrophobicity_profile
hydrophobicity_table
```

## Substitutions

```@docs
compare_proteins
```

## Codons

```@docs
codon_usage
compare_codons
```

## Grouping

```@docs
cluster_asvs
```

## Conservation

```@docs
conservation_analysis
henikoff_weights
```

## Gene length

```@docs
modal_length
```

## Anomaly scores

```@docs
robust_anomaly_scores
```

## ND2–COI concordance

```@docs
gene_concordance
spearman
```

## Simulated data

```@docs
simulate_asv_dataset
```
