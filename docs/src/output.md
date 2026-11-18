# Output tables and flags

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
