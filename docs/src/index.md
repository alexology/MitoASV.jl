# MitoASV.jl

Reference-free biological plausibility assessment of mitochondrial
protein-coding amplicon sequence variants (ASVs), designed for insect
metabarcoding and long-range ND2–COI amplicons.

This package is part of an experimental effort to use AI to port
bioinformatics and quantitative ecology packages to Julia, with the goal of
creating a unified environment for high performance scientific computing.

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

See also [examples/](https://github.com/alexology/MitoASV.jl/blob/main/examples/).

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
