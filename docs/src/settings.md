# Settings and complexity

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
