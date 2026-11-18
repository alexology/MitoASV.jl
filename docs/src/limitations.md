# Limitations and testing

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
