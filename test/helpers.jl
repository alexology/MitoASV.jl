# Test helpers. All sequences built here are SIMULATED, not real data.

# Small deterministic generator (independent of Julia's RNG implementation).
mutable struct LCG
    state::UInt64
end
next!(g::LCG) = (g.state = g.state * 6364136223846793005 + 1442695040888963407; g.state >> 33)
pick(g::LCG, v) = v[Int(next!(g) % length(v)) + 1]

# Sense codons of NCBI table 5 (no TAA/TAG), biased to hydrophobic residues
# like a membrane protein.
const SENSE5 = [c for c in (a * b * d for a in "TCAG" for b in "TCAG" for d in "TCAG")
                if !(c in ("TAA", "TAG"))]
const HYDRO_CODONS = ["TTA", "CTT", "ATT", "ATA", "GTT", "TTT", "GCT", "ATG", "TGT", "GGA"]

"""SIMULATED open reading frame of `n` codons (table 5)."""
function simulated_orf(n; seed = 1, hydrophobic = 0.5)
    g = LCG(UInt64(seed))
    io = IOBuffer()
    for _ in 1:n
        print(io, (next!(g) % 100) < hydrophobic * 100 ? pick(g, HYDRO_CODONS) : pick(g, SENSE5))
    end
    String(take!(io))
end

"""Replace codon `k` (1-based) of `s` by `codon`."""
setcodon(s, k, codon) = s[1:3(k - 1)] * codon * s[3k + 1:end]
setbase(s, i, b) = s[1:i - 1] * b * s[i + 1:end]

# One codon (table 5) for each amino acid, used to build SIMULATED variants.
const CODON_FOR = Dict('A' => "GCT", 'R' => "CGA", 'N' => "AAT", 'D' => "GAT", 'C' => "TGT",
                       'Q' => "CAA", 'E' => "GAA", 'G' => "GGA", 'H' => "CAT", 'I' => "ATT",
                       'L' => "TTA", 'K' => "AAA", 'M' => "ATA", 'F' => "TTT", 'P' => "CCT",
                       'S' => "TCT", 'T' => "ACT", 'W' => "TGA", 'Y' => "TAT", 'V' => "GTT")

"""Change the third base of codon `k` to a synonymous alternative (table 5), or return `nothing`."""
function try_synonymous(s, k)
    cod = s[3k - 2:3k]
    aa = translate_codon(cod)
    for b in "TCAG"
        alt = cod[1:2] * b
        alt != cod && translate_codon(alt) == aa && return setcodon(s, k, alt)
    end
    nothing
end
synonymous_change(s, k) = something(try_synonymous(s, k))

"""First codon index ≥ `from` whose third position has a synonymous alternative."""
syn_site(s, from) = findfirst(k -> k >= from && try_synonymous(s, k) !== nothing, 1:(length(s) ÷ 3))

"""Replace codon `k` by a codon for amino acid `aa`."""
aa_change(s, k, aa) = setcodon(s, k, CODON_FOR[aa])
