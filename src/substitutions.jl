# Amino acid substitution analysis between aligned proteins.

"""
    compare_proteins(seq1, seq2; radical = RadicalCriteria(), min_identity = 0.5) -> NamedTuple

Compare two **aligned** protein sequences position by position.

Positions are compared only if both sequences have the same length (an
alignment, or ungapped sequences known to be colinear, such as amplicons
from the same primers read in the same frame). Sequences of different
length raise an `ArgumentError`: align them first. Homology cannot be
proven without a reference; when identity over compared positions is
below `min_identity`, `homology_warning` is `true` and the metrics should
not be interpreted as substitutions between homologous residues.

Position handling:
- both standard residues: compared (`n_compared`); differing pairs are
  substitutions, scored with Grantham, BLOSUM62 and Kyte-Doolittle change
  (`seq2` minus `seq1`);
- a gap (`-`) against a residue: `n_gap_positions`; gap against gap is ignored;
- a stop (`*`) against a standard residue: `n_stop_differences`;
- `X`, `B`, `Z`, `J`, other symbols: `n_uncomparable`.

Returned metrics are `missing` when no substitution exists to summarise
(for example `grantham_mean` with zero substitutions); counts are integers.
`substitutions` lists each substitution.
"""
function compare_proteins(seq1::AbstractString, seq2::AbstractString;
                          radical::RadicalCriteria = RadicalCriteria(), min_identity::Real = 0.5)
    a = collect(uppercase(seq1)); b = collect(uppercase(seq2))
    length(a) == length(b) || throw(ArgumentError(
        "sequences differ in length ($(length(a)) vs $(length(b))); positions can only be " *
        "compared between aligned sequences. Align them first."))
    pos = Int[]; aa1 = Char[]; aa2 = Char[]; gr = Int[]; bl = Int[]; dkd = Float64[]; rad = Bool[]
    ncomp = 0; ngap = 0; nstop = 0; nunc = 0
    for i in eachindex(a)
        x, y = a[i], b[i]
        if x == '-' && y == '-'
            continue
        elseif x == '-' || y == '-'
            ngap += 1
        elseif is_standard_aa(x) && is_standard_aa(y)
            ncomp += 1
            if x != y
                push!(pos, i); push!(aa1, x); push!(aa2, y)
                push!(gr, grantham_distance(x, y)); push!(bl, blosum62_score(x, y))
                push!(dkd, kd_value(y) - kd_value(x)); push!(rad, is_radical(x, y, radical))
            end
        elseif (x == '*' && is_standard_aa(y)) || (y == '*' && is_standard_aa(x))
            nstop += 1
        elseif !(x == '*' && y == '*')
            nunc += 1
        end
    end
    nsub = length(pos)
    identity = ncomp == 0 ? missing : 1 - nsub / ncomp
    subs = DataFrame(position = pos, aa1 = aa1, aa2 = aa2, grantham = gr, blosum62 = bl,
                     delta_kd = dkd, radical = rad)
    (n_compared = ncomp, n_substitutions = nsub,
     prop_different = ncomp == 0 ? missing : nsub / ncomp,
     grantham_sum = nsub == 0 ? (ncomp == 0 ? missing : 0) : sum(gr),
     grantham_mean = nsub == 0 ? missing : mean(gr),
     grantham_max = nsub == 0 ? missing : maximum(gr),
     blosum62_mean = nsub == 0 ? missing : mean(bl),
     blosum62_min = nsub == 0 ? missing : minimum(bl),
     delta_kd_sum = ncomp == 0 ? missing : sum(dkd; init = 0.0),
     abs_delta_kd_sum = ncomp == 0 ? missing : sum(abs, dkd; init = 0.0),
     abs_delta_kd_mean = nsub == 0 ? missing : mean(abs, dkd),
     n_radical = nsub == 0 ? (ncomp == 0 ? missing : 0) : count(rad),
     n_gap_positions = ngap, n_stop_differences = nstop, n_uncomparable = nunc,
     identity = identity,
     homology_warning = identity === missing ? true : identity < min_identity,
     substitutions = subs)
end
