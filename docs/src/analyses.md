# Analyses

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
