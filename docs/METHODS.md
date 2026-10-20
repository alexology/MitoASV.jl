# MitoASV methods

This document records the design decisions behind each analysis, so that
results can be interpreted and reported. Source files are given in
brackets.

## 1. Sequence representation [`geneticcode.jl`, `sequences.jl`]

Each input sequence is parsed into a vector of IUPAC bit masks (A = 1, C = 2,
G = 4, T/U = 8, ambiguity codes are unions, gap = 16). Symbols outside IUPAC
are read as `N` and counted. The input DataFrame is never modified. Gene
segments are extracted with 1-based inclusive coordinates and reverse
complemented for `strand = :-`.

## 2. Translation and coding validation [`validation.jl`]

- NCBI tables 1, 2, 4, 5 and 9 are built from the NCBI 64-character
  strings. For every combination of three masks, the translation is the
  single amino acid shared by all resolutions, otherwise `X`.
  `n_possible_stops` counts ambiguous codons with at least one stop
  resolution.
- `frame` f means that the first complete codon starts at position f of
  the oriented segment. Bases before it (`leading_bases`) and an
  incomplete final codon (`trailing_bases`) are not translated and are
  not evidence of a frameshift.
- A stop in the last non-gap codon, with no base after it, is
  `terminal_stop` (if `terminal_stop_allowed`); any other stop is internal.
- Unknown frame: a frame is inferred only by dataset consensus, among
  unique sequences with exactly one open frame
  (`frame_min_support = 3`, `frame_consensus_min = 0.75`). All sequences
  are then read in that frame. The fewest-stop frame of an individual
  sequence is never selected, because frameshifted ASVs would otherwise
  look valid.
- Partial-gap codons in aligned input are flagged `:invalid_frame`.
  Length differences between sequences are handled by the length test (§3b).

## 3. Hydrophobicity [`hydrophobicity.jl`]

Kyte–Doolittle values are averaged over standard residues. The window
profile uses consecutive residues. In the pipeline, stops, `X` and gaps are
removed first, so that a defect already reported by the coding test does
not also change the window statistics. The 1.6 / 19-residue window
threshold is the transmembrane heuristic of Kyte & Doolittle (1982) and is
used only descriptively. Hydrophobicity enters the flags only through
within-group robust scores.

## 3b. Gene length [`length.jl`]

Usually each gene comes in its own column, and its length is simply the
number of non-gap bases. The expected length is not given by the user. It
is inferred from the data as the modal length of the **other** unique
sequences of the same gene in the same species (without `species_col`: in
the same comparison group). At least `length_min_comparators` (3) are
required, and the mode must be unique; otherwise the test is not
evaluated. A difference larger than `length_tolerance` (0 nt) raises
`:length_outlier`.

The position of the difference is checked against the closest sequence of
the expected length (at most 50 are examined):

- `terminal`: anchoring at the 5′ or 3′ end gives a mismatch fraction ≤
  `cluster_threshold`;
- `internal`: no end-anchored match, but the edit distance is ≤
  threshold × length + |difference|;
- `undetermined`: no sequence close enough to decide.

Only an internal difference that is not a multiple of 3 also raises
`:invalid_frame`. A terminal difference, such as an incomplete last codon
or a trimmed end, is not evidence of a frameshift. An indel very close to
an end can look terminal.

## 4. Grouping [`clustering.jl`, `distances.jl`]

- Distances:
  - equal lengths: p-distance with pairwise deletion of gaps and
    ambiguous bases;
  - unequal lengths: edit distance of the ungapped sequences divided by
    the longer length. It uses a banded Ukkonen algorithm whose band
    starts at |Δlength| + 8 and doubles up to the cutoff.
  - Protein distances: p-distance over standard residues, equal lengths
    only.
- Connected components: the edges at or below the threshold are computed
  in parallel and merged in a fixed order; no matrix is stored. Chaining
  is possible and is documented in the output (`nn_distance`).
- Hierarchical clustering: nearest-neighbour chain with Lance–Williams
  updates (single, complete or average linkage), cut at the threshold.
  Distances above `max(0.5, 2 × threshold)` are stored at that cap.
- Identical sequences always share a unit; `missing` sequences get no
  group.

## 4a. Groups per gene

With one column per gene and no whole-amplicon column (or
`cluster_on = :per_gene`), each gene is clustered on its own sequences.
Groups are therefore gene-specific, and a row lacking one gene is still
grouped for the others. With a whole-amplicon `sequence_col`, all genes
share groups built on the amplicon. The `clusters` table has one row per
input row and gene.

## 4b. Species [`pipeline_groups.jl`, `concordance.jl`]

With `species_col`, units are distinct (species, sequence) pairs and
clustering runs separately within each species. Comparison groups
therefore never contain two species, and neither do the within-group
statistics (§5–6). With `clustering = false` each species is one group;
with `cluster_col` the group is the (species, user group) pair. Rows with
a missing species have no group (`:missing_species`). The ND2–COI
concordance (§7) uses only ASVs of the same species
(`concordance_scope = :species`). `:all` restores dataset-wide comparisons
and is needed to find interspecific chimeras, whose second parent belongs
to another species. The reading-frame consensus (§2) stays dataset-wide
because the frame is set by the primers.

## 5. Within-group comparisons [`pipeline_groups.jl`, `conservation.jl`, `codons.jl`]

For each gene and group, the members are the **unique** gene sequences
with a resolved frame. Only members of the modal length (equal number of
columns) are compared position by position; the others are reported as
`:not_comparable`. At least `min_group_size` comparable members are
required.

- **Leave-one-out consensus.** Weighted residue counts (L × 20) and codon
  counts (64 × codons) are built once per group. For each member, its own
  contribution is subtracted. A position has a consensus only if at least
  `min_group_size − 1` comparators have an informative residue there and
  the most frequent residue is not tied.
- **Substitutions** against that consensus are scored with Grantham
  (1974), BLOSUM62 and the Kyte–Doolittle change. They are called radical
  according to `RadicalCriteria`; the default is Grantham ≥ 101 (Li et al.
  1984 classes).
- **Conservation** is 1 − H/log2 20, where H is the Shannon entropy of the
  leave-one-out frequencies of the 20 residues. Conservation-weighted
  Grantham and hydropathy changes multiply each substitution by that
  conservation.
- **Radical changes at invariant sites**: the consensus residue has
  leave-one-out frequency 1.
- **Weights.** By default each unique sequence weighs 1. With `:henikoff`,
  position-based weights (Henikoff & Henikoff 1994) are scaled to mean 1.
  With `:abundance`, weights are summed read counts; reads are not
  independent observations.
- **Codon comparisons** use the leave-one-out consensus codon at each
  position:
  - differences are counted by codon position;
  - synonymous and nonsynonymous counts follow Nei & Gojobori (1986);
  - codons that differ at 2–3 positions are averaged over all mutational
    pathways that avoid stop codons;
  - pairs involving a stop codon, or without a valid pathway, are
    `n_unclassified_codons`;
  - synonymous sites count changes to stop codons as nonsynonymous sites.
- **dN/dS** (optional, `estimate_dnds`) uses Jukes–Cantor-corrected pN and
  pS against the consensus. It is estimated only with ≥ 30 codons, at
  least one synonymous difference and p < 0.75; otherwise it is `missing`
  with a reason. It is never the raw ratio of difference counts.

## 6. Robust anomaly scores [`anomaly.jl`]

For a metric x in a gene × group (unique sequences):

    score = (x − median(others)) / dispersion(others)

The comparators are the other members (leave-one-out) when the group has
≤ `loo_dispersion_max` (200) members, otherwise all members. The
dispersion is 1.4826 × MAD. If the MAD is 0, it is 1.2533 × the mean
absolute deviation instead. In both cases it is never below a
**resolution floor**:

| Metric | Direction | Floor |
|---|---|---|
| `mean_kd` | two-sided | 4.5 / median number of residues (half the largest single-substitution change) |
| `frac_windows_above` | two-sided | 1 / median number of windows |
| `grantham_sum`, `cw_grantham_sum` | upper | 100 (one moderately radical substitution) |
| `diff_pos12`, `nonsyn_diffs` | upper | 1 (one difference) |
| `gc3` | two-sided | 1 / number of codons (one site) |

The floor prevents extreme scores that come only from zero dispersion. If
the dispersion is still 0, a value equal to the median scores 0 and any
other value gets `missing`. The default threshold is 3.5 (Iglewicz &
Hoaglin's modified z-score convention). It can be set per gene, and it
always applies to standardised scores, never to raw ND2/COI values.

The rule `n_radical_invariant ≥ 1` also raises `:grantham_outlier`; it is
recorded with `method = "rule"`.

## 7. ND2–COI concordance [`concordance.jl`]

Units are distinct (gene 1, gene 2) sequence pairs present in both genes,
taken within each species (`concordance_scope = :species`, with a species
column) or across the dataset (`:all`). All statistics below are computed
separately in each scope.
For each gene, nucleotide distances (unequal-length pairs: edit distance
capped at 0.25, beyond which the pair is not comparable) and protein
p-distances are computed into Float32 matrices. Only rank- or
neighbourhood-based statistics are used, because the genes evolve at
different rates:

- global Spearman correlation of all pairwise distances (nucleotide and
  protein). No p-value is given, because pairs are not independent;
- per ASV:
  - the nearest neighbour in each gene;
  - Jaccard overlap of the k nearest neighbours (ties included);
  - row-wise Spearman correlation;
- the two-parent model: single = min over j of
  (l1·d1(i,j) + l2·d2(i,j)) / (l1 + l2); two = the same with the best
  ASV taken separately for each gene. The improvement single − two is
  ≥ 0, and is large for a chimera between the two parents (the UCHIME
  idea applied to gene blocks);
- the minimum parent abundance divided by the ASV's abundance, when
  abundances are given;
- discordant pairs: |percentile rank in gene 1 − percentile rank in
  gene 2| ≥ `pair_discordance_threshold`;
- optional: UPGMA trees per gene and the Spearman correlation of their
  cophenetic distances.

`:gene_discordance` is raised when the improvement is ≥
`chimera_min_improvement` and the two parents differ. It is a possible
chimera indicator only.

## 8. Composite reporting [`reporting.jl`]

Each test is recorded as evaluated, not evaluated (with a reason) or not
applicable. Anomaly flags belong to four evidence classes: coding
integrity, biochemical, evolutionary and gene concordance. Comparative
flags on a gene that also has a coding-integrity flag are dependent: they
are listed but not counted as independent. The classification is
descriptive, and missing tests lead to "Insufficient evidence", never to
"No detected anomaly".
