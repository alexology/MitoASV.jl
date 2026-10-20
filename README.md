# MitoASV.jl

Reference-free biological plausibility assessment of mitochondrial
protein-coding amplicon sequence variants (ASVs), designed for insect
metabarcoding and long-range ND2–COI amplicons.

MitoASV looks for ASVs that may be artifacts using only the ASVs
themselves. It checks the reading frame, hydrophobicity, amino acid
substitutions, conservation within groups of related ASVs, codon-level
patterns and ND2–COI concordance. It needs no reference sequences, no
taxonomic assignment and no database.

> **Plausibility is not authenticity.** MitoASV reports interpretable
> scores, diagnostic flags and the statistics behind them. It never
> classifies an ASV as genuine or artifactual, and an ASV is never flagged
> only for being unusual. Without external references, read-level data or
> independent experimental evidence, no method can definitively separate
> genuine mitochondrial variants from sequencing errors, NUMTs and
> chimeras. MitoASV cannot either.

## Contents

- [Installation](#installation)
- [Minimal example](#minimal-example)
- [Comparisons within species](#comparisons-within-species)
- [Single-gene analysis](#single-gene-analysis)
- [Combined ND2–COI analysis](#combined-nd2coi-analysis)
- [What is computed](#what-is-computed)
- [Output tables and columns](#output-tables-and-columns)
- [Flags and their interpretation](#flags-and-their-interpretation)
- [Settings](#settings)
- [Computational complexity](#computational-complexity)
- [Scientific limitations](#scientific-limitations)
- [Testing and simulated data](#testing-and-simulated-data)
- [References](#references)

Method details are in [docs/METHODS.md](docs/METHODS.md).

## Installation

MitoASV requires Julia 1.10 or later. It depends only on `DataFrames` and
the `Statistics` standard library. The package is not registered; install
it from its folder or repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/alexology/MitoASV.jl")
# or, for a local copy: Pkg.develop(path = "path/to/MitoASV.jl")
```

Run the tests with `julia --project=. -e 'using Pkg; Pkg.test()'` from the
package folder (about one minute).

## Minimal example

```julia
using DataFrames, MitoASV

# One column per gene, chosen by the user, plus the species of each ASV
df = DataFrame(
    ASV     = ["ASV1", "ASV2", "ASV3"],
    species = ["Sp. a", "Sp. a", "Sp. a"],
    nd2     = ["ATT...", "ATT...", "ATT..."],     # your ND2 sequences
    coi     = ["AAC...", "AAC...", "AAC..."],     # your COI sequences
)

results = assess_asvs(df; sequence_col = nothing, id_col = :ASV,
                      genes = Dict(:ND2 => :nd2, :COI => :coi),   # gene => column
                      frame = 1, genetic_code = 5, species_col = :species)

results.summary          # one descriptive classification per input row
results.lengths          # each gene's length against the expected length
results.flags            # every flag, with the metric and threshold behind it

result_df = assess_asvs_df(df; sequence_col = nothing, id_col = :ASV,
                           genes = Dict(:ND2 => :nd2, :COI => :coi), frame = 1,
                           species_col = :species)                # one wide row per ASV
```

A single gene works the same way (`genes = [:COI => :coi]`), or with
`sequence_col = :coi, gene = :COI`.

A complete, runnable example on a simulated dataset:

```julia
using MitoASV
sim = simulate_asv_dataset(seed = 1)      # SIMULATED ND2–COI amplicons, not real data
res = assess_asvs(sim.data; sequence_col = nothing, id_col = :ASV, genes = sim.gene_columns,
                  frame = 1, species_col = :species, abundance_col = :abundance,
                  sample_col = :sample)
res                                        # prints counts per classification
```

See also [examples/](examples/).

### Input rules

- **One column per gene** is the expected input: `genes = Dict(:ND2 => :nd2,
  :COI => :coi)` (gene name => column; a vector of pairs keeps the order).
  A gene can also be `(column = :coi, frame = 1, strand = :+)`; `frame` and
  `strand` given as arguments apply to genes without their own. Rows may
  lack some genes (`missing`): each gene is analysed on its own.
- Any column names: `id_col`, `species_col`, `abundance_col`, `sample_col`,
  `replicate_col`, `cluster_col`. `sequence_col` is only needed for a
  whole-amplicon column (genes given by coordinates); set it to `nothing`
  with gene columns.
- Sequences can be strings or BioSequences objects (`LongDNA`, `LongRNA`;
  anything whose `string` is a nucleotide string). IUPAC ambiguity codes,
  gaps (`-`, `.`) and `missing` are allowed; `U` is read as `T`.
- **The input is never modified and no row is removed.** Symbols outside
  IUPAC are counted (`n_invalid_symbols`, flag `:invalid_symbols`) and
  read as `N` in the analysis copy only. Missing and empty sequences get
  flags and `missing` metrics. Every output table has the input row index
  (`row`) and identifier (`id`).
- Long format is supported: the same sequence may appear in several rows
  (samples, replicates). Each distinct sequence is analysed once, and the
  results are copied to all of its rows.

## Comparisons within species

Give the species of each ASV (for example the `species` column assigned by
a taxonomic classifier) with `species_col`. **Every comparison then stays
within a species**, because comparing variants of different species as if
they were alleles of one gene makes no sense:

```julia
res = assess_asvs(df; sequence_col = :sequence, id_col = :ASV,
                  species_col = :species, genes = genes)
```

- Clustering runs separately in each species, so a comparison group never
  contains two species. Within a species it still separates deep lineages,
  NUMT-like clusters and misassigned sequences (`cluster_threshold`).
  With `clustering = false`, each species is one comparison group. With
  `cluster_col`, the group is the (species, user group) pair.
- Consensus, conservation, codon comparisons and robust anomaly scores are
  computed within these groups, so they are always within species.
- The expected length of each gene is the modal length of the other
  sequences of that gene **in the same species** (`results.lengths`).
- The ND2–COI concordance compares an ASV only with ASVs of the same
  species (`concordance_scope = :species`, the default). An
  **interspecific** chimera (ND2 from one species, COI from another) has
  its second parent in another species. Within species it gets
  `Insufficient evidence` (it usually has no comparators), never `No
  detected anomaly`. To look for such chimeras, add a run with
  `concordance_scope = :all`.
- A sequence assigned to two species is analysed separately in each.
- Rows with a missing species get the `:missing_species` flag and no
  comparison group, so their comparative tests are not evaluated.
- Not affected by species: the per-sequence checks (frame, stops,
  hydrophobicity values) and the reading-frame consensus, which depends
  on the primers, not on the species.

Without `species_col`, groups come from clustering the whole dataset. A
fixed threshold can then merge closely related species or split one, so
giving the species is recommended whenever it is known.

## Single-gene analysis

```julia
res = assess_asvs(df; sequence_col = :sequence, gene = :COI,
                  frame = 1,          # 1, 2 or 3; `nothing` = unknown
                  strand = :+,        # :- if the sequences are reverse complemented
                  genetic_code = 5)   # NCBI table 5, invertebrate mitochondrial
```

**Give the reading frame if you know it.** With `frame = nothing`, all
three frames are evaluated. A frame is inferred only from a *dataset
consensus*: at least 3 unique sequences with exactly one open frame, at
least 75 % of them in the same frame. Every sequence is then read in that
frame. A sequence is never moved to its own fewest-stop frame, so a
frameshifted ASV gets flagged instead of being rescued. Without a
consensus, frame-dependent results are `missing` and the
`:frame_unresolved` flag explains why. Sequences need not start with ATG
or end with a stop codon. Bases before the frame and an incomplete final
codon are reported (`leading_bases`, `trailing_bases`) and are never
treated as frameshifts.

## Combined ND2–COI analysis

With one column per gene (recommended):

```julia
res = assess_asvs(df; sequence_col = nothing, id_col = :ASV,
                  genes = [:ND2 => :nd2, :COI => :coi], frame = 1,
                  species_col = :species, genetic_code = 5, abundance_col = :reads)
```

Each gene is grouped on its own sequences (within species), so an ASV
missing one gene is still assessed for the other; the concordance needs
both. Alternatively, genes can be cut from a whole-amplicon column by
coordinates; the groups are then built on the whole amplicon:

```julia
res = assess_asvs(df; sequence_col = :sequence,
    genes = Dict(
        :ND2 => (start = 1,    stop = 900,  frame = 1),
        :COI => (start = 1200, stop = 2700, frame = 1),
    ),
    genetic_code = 5, abundance_col = :reads)
```

- Coordinates are **1-based and inclusive** on the input sequence as given
  (alignment columns if the input is aligned). The coordinates above are
  examples only: use those of your amplicon.
- Each gene can have its own `frame` and `strand` (`:-` means the gene is
  the reverse complement of `start:stop`).
- `cluster_on = :per_gene` forces per-gene groups, `cluster_on = :COI`
  groups all genes on COI.
- A sequence that ends inside a gene is analysed as a partial gene
  (`:gene_truncated`). If it ends before the gene starts, it gets
  `:gene_absent`.
- Every analysis is done separately for each gene. ND2 and COI are never
  compared with the same raw thresholds: anomaly scores are standardised
  within each gene and group, and `anomaly_threshold` can be set per gene,
  e.g. `Dict(:ND2 => 4.0, :COI => 3.5)`.
- The two genes are compared in `results.gene_concordance` and
  `results.discordant_pairs` (see below).

## What is computed

| Step | Method | Main output |
|---|---|---|
| Coding validation | Translation with NCBI table 5 (also 1, 2, 4, 9). Ambiguous codons are expanded: they get an amino acid only if every resolution agrees, otherwise `X`, with possible stops counted. The step also reports internal and terminal stops, partial-gap codons, GC overall and at codon positions 1–3, and amino acid composition. | `validation`, `composition` |
| Hydrophobicity | Kyte–Doolittle residue statistics and a sliding window (default 19 residues). | `hydrophobicity` |
| Gene length | Length of each gene sequence against the **expected length inferred from the data**: the modal length of the other unique sequences of that gene in the same species (or group). A difference is located at an end (terminal) or inside the gene (internal indel). | `lengths` |
| Grouping | Connected components (default) or hierarchical clustering of pairwise distances, on nucleotides (default) or proteins; groups can also be given by the user. | `clusters` |
| Substitutions | Comparison with the group's **leave-one-out** consensus: Grantham, BLOSUM62, hydropathy change, radical substitutions. | `amino_acid_changes` |
| Conservation | Leave-one-out position frequencies, Shannon entropy and conservation; unseen and rare residues; conservation-weighted Grantham and hydropathy changes; radical changes at sites invariant among the comparators. | `conservation` |
| Codons | Differences at codon positions 1, 2 and 3. Synonymous/nonsynonymous counts follow Nei–Gojobori (1986), with codons that differ at several positions averaged over all stop-free pathways. dN/dS is optional and computed only when estimable. Codon usage and RSCU are also reported. | `codon_statistics`, `codon_usage`, `group_statistics` |
| Anomalies | Robust, leave-one-out within-group scores with explicit handling of zero dispersion. | `anomalies` |
| ND2–COI concordance | Rank correlation of the gene distance matrices, nearest-neighbour and k-NN consistency per ASV, two-parent versus one-parent fit, and discordant pairs. | `gene_concordance`, `discordant_pairs` |
| Report | Flags in long format, a record of every test and its status, and a descriptive classification. | `flags`, `tests`, `summary` |

Identical sequences are collapsed before every group statistic, so
duplicates cannot inflate conservation. Groups with fewer than
`min_group_size` (default 5) unique, equal-length sequences give
`:insufficient_comparators` and no scores.

## Output tables and columns

`assess_asvs` returns an `AssessmentResult`. All its fields are
DataFrames. Per-ASV tables have one row per input row, or one row per
input row and gene, keyed by `row` and `id`. A `reason` column says why a
value is `missing`.

### `summary` (one row per input row)

| Column | Meaning |
|---|---|
| `species` | the `species_col` value (only when given) |
| `classification` | `No detected anomaly`, `One or more anomalies`, `Multiple independent anomalies` or `Insufficient evidence` (see below) |
| `n_anomaly_flags`, `anomaly_flags` | anomaly flags raised, as `GENE:flag` |
| `n_evidence_classes`, `evidence_classes` | distinct independent evidence classes among them (the composite summary) |
| `dependent_flags` | comparative flags on a gene that also has a coding-integrity defect; they are consequences of that defect and are not counted as independent evidence |
| `n_tests_evaluated`, `n_tests_not_evaluated`, `tests_not_evaluated` | test bookkeeping (details in `tests`) |
| `other_flags` | data and status flags (missing sequence, truncated gene, insufficient comparators, ...) |
| `n_rows_same_sequence`, `total_abundance`, `n_samples_detected`, `n_replicates_detected` | detection support, descriptive only and never used as evidence |

### `validation` (row × gene)

`status` (`ok`, `truncated`, `absent`, `missing_sequence`,
`empty_sequence`); `frame_source` (`user`, `consensus`, `unresolved`);
`frame`; `length_nt` (non-gap bases); `n_columns` (including gaps);
`n_ambiguous_nt`; `gc`, `gc1`, `gc2`, `gc3` (GC fraction overall and at
codon positions 1–3, ambiguous bases excluded); `leading_bases`;
`n_codons` (complete codons); `trailing_bases` (incomplete final codon,
0–2); `n_internal_stops`; `internal_stop_codons` (codon positions);
`terminal_stop`; `n_ambiguous_codons`; `n_unresolved_codons` (translated
`X`); `n_possible_stops` (ambiguous codons that may be stops);
`n_partial_gap_codons`; `n_gap_codons`; `open_frames` (open forward
frames); `three_frame_status` (`unique_open_frame`,
`multiple_open_frames`, `no_open_frame`); `n_invalid_symbols`;
`protein`.

### `composition` (row × gene)

Columns `A` … `Y`: fraction of each standard amino acid among standard
residues.

### `hydrophobicity` (row × gene)

`n_residues_scored`; `mean_kd`, `median_kd`, `sd_kd` (Kyte–Doolittle);
`frac_hydrophobic` (residues with KD > `hydrophobic_cutoff`, default 0);
`n_windows` (complete windows); `frac_windows_above` (window mean >
`window_threshold`, default 1.6); `max_window_kd`, `min_window_kd`. In
`assess_asvs`, stops, `X` and gaps are skipped, so a stop already reported
by the coding test does not distort the windows. The stand-alone
`hydrophobicity_profile` instead returns `missing` for windows that span
such residues.

### `lengths` (row × gene)

`scope` (`species`, or `group` without species); `length_nt` (non-gap
bases of the gene sequence); `modal_length_nt` (expected length: modal
length of the **other** unique sequences of the gene in the same species
or group; `missing` with fewer than `length_min_comparators` = 3 of them
or a tie); `length_difference`; `multiple_of_3`; `indel_location`
(`terminal`, `internal`, `undetermined`, for differing lengths);
`n_comparators`; `status`, `reason`.

### `amino_acid_changes` (row × gene, against the leave-one-out group consensus)

`cluster`; `n_compared` (positions with a consensus and a standard
residue); `n_substitutions`; `prop_different`; `grantham_sum`,
`grantham_mean`, `grantham_max`; `blosum62_mean`, `blosum62_min`;
`delta_kd_sum` (signed), `abs_delta_kd_sum`; `n_radical`;
`n_stop_vs_consensus`; `n_gap_vs_consensus`;
`n_positions_without_consensus` (tie, or too few informative comparators).

### `conservation` (row × gene)

`n_scored_positions`; `mean_conservation` (1 − entropy/log2 20, leave-one-out);
`n_unseen_residues` (residue absent from all comparators at that position);
`n_rare_residues` (frequency < `rare_threshold`); `cw_grantham_sum`,
`cw_grantham_max`, `cw_abs_delta_kd_sum` (each substitution weighted by
the site's leave-one-out conservation); `n_radical_invariant` (radical
substitutions at sites invariant among the comparators); `n_comparators`;
`n_effective_group` (Kish effective number from the weights).

### `codon_statistics` (row × gene, against the leave-one-out consensus codons)

`n_codons_compared`; `diff_pos1`, `diff_pos2`, `diff_pos3`, `diff_pos12`
(nucleotide differences by codon position); `n_diff_codons`;
`n_multi_diff_codons`; `syn_diffs`, `nonsyn_diffs` (Nei–Gojobori,
fractional for multi-difference codons); `prop_nonsyn`;
`n_unclassified_codons` (stop involved or no valid pathway);
`syn_sites`, `nonsyn_sites`; `pN`, `pS`, `dN`, `dS`, `dN_dS` and
`dnds_reason` (only with `estimate_dnds = true`, and only when estimable).

### `codon_usage` (gene × codon)

`codon`, `amino_acid`, `count` (over unique gene sequences),
`n_synonymous_codons`, `rscu` (`missing` for stops, single-codon amino
acids and amino acids never observed).

### `clusters` (one row per input row and gene)

`gene` (with shared groups the rows repeat for each gene); `species`; `unit` (distinct (species, sequence) pair); `cluster`;
`cluster_label` (the group label: cluster number, user label, or
(species, label) with `species_col`); `n_unique_in_cluster`; `n_rows_in_cluster`; `nn_distance`
(nearest other sequence, if within the threshold); `cluster_source`
(`connected_components`, `hierarchical`, `user`, `single_group`);
`cluster_basis`; `threshold`; `reason` (why a row has no group).

### `group_statistics` (gene × group)

`gene`; `species`; `cluster_label`; `cluster`; `n_members`, `n_comparable` (unique, equal-length); `status`, `reason`;
`modal_length_nt`; `n_codons`; `n_effective`; `var_frac_pos1..3`
(fraction of codons variable at each position); `aa_variable_frac`;
`mean_conservation`.

### `gene_concordance` (one row per input row; multi-gene runs)

`gene_pair`; `scope` (`species` or `all`); `species`; `nn_gene1_id`,
`nn_gene2_id` (nearest neighbour in each gene, within the scope); `distance_nn_gene1`, `distance_nn_gene2`; `knn_jaccard` (overlap of
the k nearest neighbours in the two genes, ties included);
`row_spearman_nt`, `row_spearman_aa` (rank correlation of this ASV's
distances to all others in the two genes); `single_parent_distance`,
`two_parent_distance`, `two_parent_improvement`;
`parents_abundance_ratio` (minimum parent abundance divided by this ASV's
abundance); `n_comparators`; `discordant` (the `:gene_discordance`
rule); `global_spearman_nt`, `global_spearman_aa` (over the scope: the
species, or the dataset); `tree_spearman`
(optional UPGMA cophenetic correlation); `status`, `reason`. In
single-gene runs the table only holds `status = not_applicable`.

### `discordant_pairs`

ASV pairs whose distance percentiles differ by at least
`pair_discordance_threshold` between the genes, within the scope. Columns:
`species`, `row_i`, `id_i`,
`row_j`, `id_j`, `distance_gene1`, `distance_gene2`, `percentile_gene1`,
`percentile_gene2`, `discordance`.

### `anomalies` (row × gene × metric)

`metric`; `test`; `method` (`robust_score_loo`, `robust_score` or `rule`);
`direction` (`two_sided` or `upper`); `value` (raw value);
`group_median`; `group_dispersion`; `dispersion_method` (`mad`,
`mean_absolute_deviation`, `resolution_floor`, `zero_dispersion`);
`robust_score`; `threshold`; `is_outlier`; `reliable` (whether a score
could be estimated); `n_comparators`; `reason`.

### `flags`, `tests`, `parameters`

- `flags`: `gene`, `flag`, `flag_type` (`anomaly`, `data`, `status`),
  `evidence_class`, `metric`, `value`, `threshold` and a human-readable
  `detail`.
- `tests`: for every row, gene and test (`coding`, `hydrophobicity`,
  `grantham`, `codon_position`, `nonsynonymous`, `length`, `gene_concordance`), the
  `status` (`evaluated`, `not_evaluated`, `not_applicable`) and a `reason`.
- `parameters`: every setting used, the frames used and how they were
  obtained, the column names, the number of threads and the package
  version.

### `assess_asvs_df`

Returns one row per input row with:

- the original columns;
- the summary columns;
- `<GENE>_cluster` for each gene;
- `flag_<name>` columns;
- the main metrics of each gene, prefixed with the gene name (for example
  `COI_mean_kd`, `ND2_n_internal_stops`);
- `concordance_*` columns.

The flag columns are `true` (raised), `false` (test evaluated for every
applicable gene and not raised) or `missing` (not evaluated). Names that
clash with input columns get a `mito_` prefix.

## Flags and their interpretation

| Flag | Class | Raised when | Interpretation |
|---|---|---|---|
| `:internal_stop` | coding integrity | a stop codon before the last codon in the frame used | strong sign of a sequencing error, frameshift or pseudogene (NUMT). Check the frame first |
| `:invalid_frame` | coding integrity | no complete codon; a codon with 1–2 gap columns; or an **internal** indel whose length is not a multiple of 3 (length test) | frameshift-like indel (or misalignment). A length difference at an end, such as an incomplete terminal codon, is never treated as a frameshift |
| `:length_outlier` | coding integrity | the gene length differs from the expected length (modal length of the other sequences of the species) by more than `length_tolerance` (0 nt) | in-frame indel, truncation, extension or frameshift (see `indel_location`, `multiple_of_3`); real codon indels exist but are rare within a species |
| `:hydrophobicity_outlier` | biochemical | \|robust score\| of `mean_kd` or `frac_windows_above` ≥ threshold | weak evidence; hydrophobicity does not prove function |
| `:grantham_outlier` | biochemical | robust score of `grantham_sum` or `cw_grantham_sum` ≥ threshold, **or** ≥ 1 radical substitution at a site invariant among the comparators | unusual physico-chemical change compared with related ASVs |
| `:codon_position_outlier` | evolutionary | robust score of `diff_pos12` (upper) or `gc3` (two-sided) ≥ threshold | differences concentrated at codon positions 1–2, or unusual GC3; random errors fall at all positions, while real variation is mostly at third positions |
| `:nonsynonymous_outlier` | evolutionary | robust score of `nonsyn_diffs` ≥ threshold | more amino acid-changing differences than the group |
| `:gene_discordance` | gene concordance | the two genes have different nearest neighbours **and** two parents fit at least `chimera_min_improvement` better than one | possible chimera indicator, not evidence |
| `:insufficient_comparators` | status | fewer than `min_group_size` unique comparable sequences | no comparative test possible: **not** a pass |
| `:frame_unresolved`, `:not_comparable` | status | no frame consensus; length differs from the group's modal length | the tests concerned are not evaluated |
| `:missing_sequence`, `:missing_species`, `:empty_sequence`, `:invalid_symbols`, `:gene_absent`, `:gene_truncated`, `:possible_stop` | data | as named | describe the input |

### Classification (`summary.classification`)

- **Multiple independent anomalies**: anomaly flags from at least two
  evidence classes. Comparative flags on a gene with a coding-integrity
  defect do not count as independent.
- **One or more anomalies**: at least one anomaly flag.
- **No detected anomaly**: no anomaly flag **and** every applicable test
  was evaluated.
- **Insufficient evidence**: no anomaly flag, but at least one applicable
  test could not be evaluated (small group, unresolved frame, missing
  gene, disabled analysis, ...). **A missing test result is never
  counted as passed.**

The composite summary is the number of independent evidence classes. It
is *not* a probability that an ASV is an artifact, and heterogeneous raw
metrics are never summed. All individual evidence stays accessible in
`flags`, `anomalies` and the per-analysis tables.

## Settings

All settings are keyword arguments of `assess_asvs` and are recorded in
`results.parameters`; see `?assess_asvs` for the full list. The main ones:

| Setting | Default | Note |
|---|---|---|
| `genetic_code` | 5 | NCBI table |
| `frame`, `strand` | `nothing`, `:+` | per gene in `genes` |
| `terminal_stop_allowed` | `true` | `false` counts a final stop as internal (amplicons inside a gene) |
| `hydrophobicity_window`, `window_threshold` | 19, 1.6 | descriptive, not validity thresholds |
| `cluster_method` | `:connected_components` | or `:hierarchical`, `:none` (all ASVs form one group) |
| `cluster_threshold` | 0.03 | an example value, **not** a species boundary |
| `cluster_level` | `:nucleotide` | or `:protein` |
| `min_group_size` | 5 | unique sequences needed for comparative tests |
| `radical` | `RadicalCriteria(grantham_min = 101)` | configurable, see `?RadicalCriteria` |
| `conservation_weights` | `nothing` | `:henikoff` or `:abundance` |
| `anomaly_threshold` | 3.5 | number or `Dict` gene ⇒ threshold |
| `estimate_dnds` | `false` | NG86 + Jukes–Cantor, only when estimable |
| `comparative`, `concordance`, `clustering` | `true` | turn off expensive steps |
| `length_tolerance`, `length_min_comparators` | 0, 3 | allowed length difference (nt); other unique sequences needed to infer the expected length |
| `species_col` | `nothing` | comparisons within species (recommended) |
| `concordance_scope` | `:species` | `:all` also compares across species (interspecific chimeras) |
| `compare_gene_trees` | `false` | UPGMA cophenetic comparison |
| `threaded` | `true` | results do not depend on the number of threads |

Stand-alone functions:

- `translate_sequence`, `coding_validation`, `frame_candidates`;
- `hydrophobicity_profile`, `hydrophobicity_metrics`;
- `compare_proteins`, `compare_codons`, `codon_usage`;
- `cluster_asvs`, `conservation_analysis`;
- `robust_anomaly_scores`, `gene_concordance`;
- `grantham_distance`, `blosum62_score`, `is_radical`.

## Computational complexity

n = distinct sequences, L = sequence length, m = members of a group.

| Step | Time | Memory |
|---|---|---|
| Parsing, translation, hydrophobicity | O(n L), each distinct sequence once | O(n L) |
| Connected components | O(n² L) worst case; the Hamming comparison stops as soon as the threshold is exceeded, and unequal lengths use a banded edit distance whose band starts narrow | edges only, no matrix |
| Hierarchical clustering | O(n² L) distances + O(n²) nearest-neighbour chain | n × n matrix (limited by `max_hierarchical = 3000`) |
| Comparative analyses | O(m L) per group, with leave-one-out by count subtraction | O(L) per group |
| Robust scores | O(m² log m) with leave-one-out (m ≤ 200), otherwise O(m log m) | O(m) |
| Concordance | O(n² L) | four n × n Float32 matrices (limited by `max_concordance_asvs = 3000`) |

Distances and per-sequence metrics use Julia threads. Results are
written to fixed slots or merged in a fixed order, so they are identical
for any number of threads (start Julia with `-t auto`). As a guide,
2,105 simulated 2.6 kb ND2–COI ASVs (1,798 distinct) took about 7 s with
4 threads (4 s without concordance) on a laptop.

## Scientific limitations

- **No reference, no proof.** Without external references, read-level
  information (quality scores, read counts per replicate, error models)
  or independent evidence (e.g. RT-PCR, morphology), MitoASV cannot
  definitively tell genuine mitochondrial variants from sequencing
  errors, **NUMTs** and **chimeras**. Recent NUMTs can be fully
  open-reading and biochemically ordinary; real variants can be rare and
  unusual.
- **Relative tests depend on groups.** Comparative scores are only as
  good as the comparison groups. With `species_col` they are only as good
  as the species assignment: a misassigned sequence is compared with the
  wrong species (clustering within the species usually isolates it, which
  gives insufficient comparators rather than a false score). Connected components can chain
  divergent sequences; any fixed threshold can split a species or merge
  several. Groups that are too small give no scores (insufficient
  comparators) rather than unreliable ones.
- **ASVs are not independent.** Related ASVs share ancestry, so
  conservation estimated from them is biased toward the sampled
  diversity. Identical sequences are collapsed and leave-one-out
  estimates are used, but closely related sequences still dominate.
  Henikoff weights reduce the effect without removing it.
- **Homology is assumed, not proven.** Positions are compared only
  between equal-length sequences read in the same frame within a group
  (or columns of an aligned input). Sequences of different length are
  not compared without an alignment.
- **Hydrophobicity is weak evidence.** Thresholds are not universal for
  ND2 or COI, and a plausible profile does not prove function.
- **Codon statistics.** Nei–Gojobori pathway averaging is an
  approximation. Between closely related sequences dN/dS reflects
  polymorphism, not selection, and is off by default.
- **Concordance is indirect.** Animal mtDNA does not normally
  recombine, but incomplete sampling, homoplasy, rate differences
  between genes, NUMTs and errors also create discordance. A discordant
  ASV is a candidate for checking (e.g. with UCHIME-type tools and read
  data), not a confirmed chimera.
- **Flags are diagnostics.** Thresholds (3.5 robust units, Grantham ≥ 101,
  ...) are conventions recorded in the output. Change them to fit your
  data and report what you used.

## Testing and simulated data

`test/runtests.jl` covers:

- translation with table 5, ambiguous nucleotides, internal stops,
  partial sequences and invalid frames;
- hydrophobicity, Grantham distances, synonymous/nonsynonymous
  classification and clustering;
- conservation, robust scores and gene concordance;
- missing values, small groups, identical sequences and zero-variance
  groups.

`simulate_asv_dataset` builds a **simulated** ND2–spacer–COI dataset
(random open reading frames under table 5, *not real sequences*) with the
following categories:

- plausible variants;
- a single conservative nonsynonymous change;
- a radical substitution;
- an internal stop;
- a 1-nt deletion;
- an ND2–COI recombinant;
- a divergent but plausible lineage;
- a second, unrelated species (`species` column: `Simulatus alpha`,
  `Simulatus beta`).

All simulated rows have `simulated = true`.

## References

- Grantham R. (1974) Amino acid difference formula to help explain protein evolution. *Science* 185:862–864.
- Henikoff S., Henikoff J.G. (1992) Amino acid substitution matrices from protein blocks. *PNAS* 89:10915–10919.
- Henikoff S., Henikoff J.G. (1994) Position-based sequence weights. *J. Mol. Biol.* 243:574–578.
- Iglewicz B., Hoaglin D.C. (1993) *How to Detect and Handle Outliers*. ASQC Quality Press.
- Kryazhimskiy S., Plotkin J.B. (2008) The population genetics of dN/dS. *PLoS Genetics* 4:e1000304.
- Kyte J., Doolittle R.F. (1982) A simple method for displaying the hydropathic character of a protein. *J. Mol. Biol.* 157:105–132.
- Li W.-H., Wu C.-I., Luo C.-C. (1984) Nonrandomness of point mutation as reflected in nucleotide substitutions in pseudogenes. *J. Mol. Evol.* 21:58–71.
- Nei M., Gojobori T. (1986) Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. *Mol. Biol. Evol.* 3:418–426.
- Sharp P.M., Li W.-H. (1987) The codon adaptation index. *Nucleic Acids Res.* 15:1281–1295.
- Ukkonen E. (1985) Algorithms for approximate string matching. *Information and Control* 64:100–118.
- Edgar R.C. et al. (2011) UCHIME improves sensitivity and speed of chimera detection. *Bioinformatics* 27:2194–2200.

## License

MIT (see `LICENSE`).
