# SIMULATED data for examples and tests. Nothing here is real sequence data.

mutable struct _SimRNG
    state::UInt64
end
_next!(g::_SimRNG) = (g.state = g.state * 0x5851f42d4c957f2d + 0x14057b7ef767814f; g.state >> 33)
_rand(g::_SimRNG, n::Integer) = Int(_next!(g) % n) + 1          # 1..n
_pick(g::_SimRNG, v) = v[_rand(g, length(v))]

const _SIM_SENSE = [codon_string(c) for c in 1:64 if genetic_code(5).table[c] != '*']
const _SIM_HYDRO = ["TTA", "CTT", "ATT", "ATA", "GTT", "TTT", "GCT", "ATG", "TGT", "GGA", "CTA", "TTC"]

function _sim_orf(g, ncod)
    io = IOBuffer()
    for _ in 1:ncod
        print(io, _rand(g, 100) <= 55 ? _pick(g, _SIM_HYDRO) : _pick(g, _SIM_SENSE))
    end
    String(take!(io))
end

_sim_setcodon(s, k, c) = s[1:3(k - 1)] * c * s[3k + 1:end]
_sim_aa(c) = translate_codon(c)

"""Synonymous third-position change at codon k (table 5), or `nothing`."""
function _sim_syn(g, s, k)
    cod = s[3k - 2:3k]; aa = _sim_aa(cod)
    alts = [cod[1:2] * b for b in "TCAG" if cod[1:2] * b != cod && _sim_aa(cod[1:2] * b) == aa]
    isempty(alts) ? nothing : _sim_setcodon(s, k, _pick(g, alts))
end

function _sim_syn_changes(g, s, n; avoid = Int[])
    ncod = length(s) ÷ 3; done = 0; tries = 0
    while done < n && tries < 100n
        tries += 1
        k = _rand(g, ncod)
        k in avoid && continue
        t = _sim_syn(g, s, k)
        t === nothing && continue
        s = t; done += 1
    end
    s
end

const _SIM_CONSERVATIVE = Dict('I' => "GTT", 'V' => "ATT", 'L' => "ATT", 'F' => "TTA",
                               'A' => "TCT", 'S' => "GCT", 'T' => "TCT", 'M' => "TTA")

"""
    simulate_asv_dataset(; seed = 1, nd2_codons = 330, coi_codons = 510, spacer = 120,
                         n_variants = 10) -> NamedTuple

Build a **SIMULATED** long-range ND2–spacer–COI amplicon dataset for
examples and tests (random open reading frames under NCBI table 5; not
real sequences). Returns `data` (a DataFrame with `ASV`, `sequence` (whole
amplicon), `species`, `category`, `lineage`, `abundance`, `sample`,
`simulated = true`, and one column per gene: `ND2`, `COI`), `genes`
(coordinates on `sequence`), `gene_columns` (`[:ND2 => :ND2, :COI => :COI]`,
for `assess_asvs(data; sequence_col = nothing, genes = gene_columns, frame = 1)`)
and `description`.

Species: `Simulatus alpha` holds lineages A and B (two deep intraspecific
mitochondrial lineages) and categories 1–6; `Simulatus beta` (category 8)
is an unrelated species. Categories:
1. `plausible_variant`: lineage A variants with 1–3 synonymous changes;
2. `single_nonsynonymous`: one conservative amino acid change in COI;
3. `radical_substitution`: one radical change at a COI site invariant in A;
4. `internal_stop`: a stop codon inside COI;
5. `frameshift`: a 1-nt deletion inside COI;
6. `recombinant`: ND2 (and spacer) from lineage A + COI from lineage B;
7. `divergent_plausible`: lineage B, ~25 % of codons changed (≈ 8 % of
   nucleotides; synonymous third positions and a few conservative changes),
   with its own variants;
8. `other_species`: variants of an independent random ORF pair, species
   `Simulatus beta`.
"""
function simulate_asv_dataset(; seed::Integer = 1, nd2_codons::Integer = 330,
                              coi_codons::Integer = 510, spacer::Integer = 120,
                              n_variants::Integer = 10)
    g = _SimRNG(UInt64(seed) * 0x9e3779b97f4a7c15 + 1)
    nd2A = _sim_orf(g, nd2_codons); coiA = _sim_orf(g, coi_codons)
    spA = String([_rand(g, 100) <= 70 ? _pick(g, ('A', 'T')) : _pick(g, ('C', 'G')) for _ in 1:spacer])
    # lineage B: divergent but plausible
    nd2B = _sim_syn_changes(g, nd2A, nd2_codons ÷ 4)
    coiB = _sim_syn_changes(g, coiA, coi_codons ÷ 4)
    for _ in 1:3
        k = _rand(g, coi_codons); a = _sim_aa(coiB[3k - 2:3k])
        haskey(_SIM_CONSERVATIVE, a) && (coiB = _sim_setcodon(coiB, k, _SIM_CONSERVATIVE[a]))
    end
    spB = _sim_syn_changes(g, spA, 0)
    rows = NamedTuple[]
    add!(id, seq, cat, lin, ab, species = "Simulatus alpha") =
        push!(rows, (ASV = id, sequence = seq, species = species, category = cat, lineage = lin,
                     abundance = ab, sample = _rand(g, 2) == 1 ? "S1" : "S2", simulated = true))
    varsA = String[]; varsB = String[]
    for i in 1:n_variants
        n2 = _sim_syn_changes(g, nd2A, _rand(g, 3) - 1)
        c = _sim_syn_changes(g, coiA, _rand(g, 3) - 1)
        s = n2 * spA * c
        push!(varsA, s)
        add!("SIM_A$(lpad(i, 2, '0'))", s, "plausible_variant", "A", 100 * _rand(g, 10))
    end
    for i in 1:max(6, n_variants ÷ 2)
        s = _sim_syn_changes(g, nd2B, _rand(g, 3) - 1) * spB * _sim_syn_changes(g, coiB, _rand(g, 3) - 1)
        push!(varsB, s)
        add!("SIM_B$(lpad(i, 2, '0'))", s, "divergent_plausible", "B", 100 * _rand(g, 10))
    end
    coi_start = 3nd2_codons + spacer + 1
    base = nd2A * spA
    # 2. conservative change
    k = findfirst(k -> haskey(_SIM_CONSERVATIVE, _sim_aa(coiA[3k - 2:3k])), 50:coi_codons) + 49
    add!("SIM_NS1", base * _sim_setcodon(coiA, k, _SIM_CONSERVATIVE[_sim_aa(coiA[3k - 2:3k])]),
         "single_nonsynonymous", "A", 40)
    # 3. radical change: a hydrophobic residue replaced by Asp (Grantham ≥ 126)
    k = findfirst(k -> _sim_aa(coiA[3k - 2:3k]) in "ILVFMCWYA", (3coi_codons ÷ 4):coi_codons) +
        3coi_codons ÷ 4 - 1
    add!("SIM_RAD1", base * _sim_setcodon(coiA, k, "GAT"),
         "radical_substitution", "A", 12)
    # 4. internal stop
    add!("SIM_STOP1", base * _sim_setcodon(coiA, coi_codons ÷ 3, "TAA"), "internal_stop", "A", 6)
    # 5. frameshift-like 1-nt deletion
    p = 3(coi_codons ÷ 2)
    add!("SIM_FS1", base * coiA[1:p - 1] * coiA[p + 1:end], "frameshift", "A", 5)
    # 6. ND2 of A + COI of B
    add!("SIM_CHIM1", varsA[1][1:(coi_start - 1)] * varsB[1][coi_start:end], "recombinant", "A×B", 8)
    # 8. a second species: independent random ORFs (unrelated to alpha)
    nd2C = _sim_orf(g, nd2_codons); coiC = _sim_orf(g, coi_codons)
    spC = String([_rand(g, 100) <= 70 ? _pick(g, ('A', 'T')) : _pick(g, ('C', 'G')) for _ in 1:spacer])
    for i in 1:max(6, n_variants ÷ 2)
        s = _sim_syn_changes(g, nd2C, _rand(g, 3) - 1) * spC * _sim_syn_changes(g, coiC, _rand(g, 3) - 1)
        add!("SIM_C$(lpad(i, 2, '0'))", s, "other_species", "C", 100 * _rand(g, 10), "Simulatus beta")
    end
    data = DataFrame(rows)
    genes = Dict(:ND2 => (start = 1, stop = 3nd2_codons, frame = 1),
                 :COI => (start = coi_start, stop = coi_start + 3coi_codons - 1, frame = 1))
    # one column per gene, as in a typical input table
    data.ND2 = [s[1:3nd2_codons] for s in data.sequence]
    data.COI = [s[coi_start:min(end, coi_start + 3coi_codons - 1)] for s in data.sequence]
    (data = data, genes = genes, gene_columns = [:ND2 => :ND2, :COI => :COI],
     description = "SIMULATED ND2–spacer–COI amplicons (NCBI table 5); random ORFs, not real sequences")
end
