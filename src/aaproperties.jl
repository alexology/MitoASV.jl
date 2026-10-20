# Amino acid property tables: Kyte-Doolittle hydropathy, Grantham distances,
# BLOSUM62, and the configurable definition of radical substitutions.

include("property_data.jl")

# Index 1..20 of a standard amino acid byte (order AMINO_ACIDS), 0 otherwise.
const AA_INDEX = let t = zeros(Int, 256)
    for (i, c) in enumerate(AMINO_ACIDS)
        t[UInt8(c) + 1] = i
    end
    t
end
@inline aa_index(c::UInt8) = @inbounds AA_INDEX[c + 1]
@inline aa_index(c::Char) = isascii(c) ? aa_index(UInt8(c)) : 0
@inline is_standard_aa(c) = aa_index(c) != 0

"""
Kyte & Doolittle (1982) hydropathy index, J. Mol. Biol. 157:105-132.
"""
const KYTE_DOOLITTLE = Dict{Char,Float64}(
    'A' => 1.8, 'R' => -4.5, 'N' => -3.5, 'D' => -3.5, 'C' => 2.5, 'Q' => -3.5,
    'E' => -3.5, 'G' => -0.4, 'H' => -3.2, 'I' => 4.5, 'L' => 3.8, 'K' => -3.9,
    'M' => 1.9, 'F' => 2.8, 'P' => -1.6, 'S' => -0.8, 'T' => -0.7, 'W' => -0.9,
    'Y' => -1.3, 'V' => 4.2)

const KD_VEC = [KYTE_DOOLITTLE[c] for c in AMINO_ACIDS]
@inline kd_value(c) = (k = aa_index(c); k == 0 ? NaN : @inbounds KD_VEC[k])

const THREE_TO_ONE = Dict(
    "Ala" => 'A', "Arg" => 'R', "Asn" => 'N', "Asp" => 'D', "Cys" => 'C', "Gln" => 'Q',
    "Glu" => 'E', "Gly" => 'G', "His" => 'H', "Ile" => 'I', "Leu" => 'L', "Lys" => 'K',
    "Met" => 'M', "Phe" => 'F', "Pro" => 'P', "Ser" => 'S', "Thr" => 'T', "Trp" => 'W',
    "Tyr" => 'Y', "Val" => 'V')

const GRANTHAM_MATRIX = let
    lines = filter(!isempty, split(GRANTHAM_CSV, '\n'))
    header = split(lines[1], ',')[2:end]
    M = fill(-1, 20, 20)
    for line in lines[2:end]
        f = split(line, ',')
        i = aa_index(THREE_TO_ONE[f[1]])
        for (h, v) in zip(header, f[2:end])
            M[i, aa_index(THREE_TO_ONE[h])] = parse(Int, v)
        end
    end
    all(>=(0), M) && M == M' || error("Grantham table is incomplete or not symmetric")
    M
end

const BLOSUM62_SYMBOLS = "ARNDCQEGHILKMFPSTWYVBZX*"
const BLOSUM62_MATRIX = let
    lines = filter(!isempty, split(BLOSUM62_TSV, '\n'))
    header = [only(h) for h in split(lines[1], '\t')]
    M = zeros(Int, 24, 24)
    for line in lines[2:end]
        f = split(line, '\t')
        i = findfirst(==(only(f[1])), BLOSUM62_SYMBOLS)
        for (h, v) in zip(header, f[2:end])
            M[i, findfirst(==(h), BLOSUM62_SYMBOLS)] = parse(Int, v)
        end
    end
    M == M' || error("BLOSUM62 table is not symmetric")
    M
end
const BLOSUM62_INDEX = let t = zeros(Int, 256)
    for (i, c) in enumerate(BLOSUM62_SYMBOLS)
        t[UInt8(c) + 1] = i
    end
    t
end

"""
    grantham_distance(a, b) -> Union{Int,Missing}

Grantham (1974) physico-chemical distance between two standard amino acids
(one-letter codes). `missing` if either is not one of the 20 standard
residues (gap, stop, `X`, `B`, `Z`, ...).
"""
function grantham_distance(a, b)
    i = aa_index(_aachar(a)); j = aa_index(_aachar(b))
    (i == 0 || j == 0) ? missing : GRANTHAM_MATRIX[i, j]
end

"""
    blosum62_score(a, b) -> Union{Int,Missing}

BLOSUM62 score. Defined for the 20 standard residues, `B`, `Z`, `X` and
`*`; `missing` for gaps and other symbols.
"""
function blosum62_score(a, b)
    ca = _aachar(a); cb = _aachar(b)
    (isascii(ca) && isascii(cb)) || return missing
    i = BLOSUM62_INDEX[UInt8(ca) + 1]; j = BLOSUM62_INDEX[UInt8(cb) + 1]
    (i == 0 || j == 0) ? missing : BLOSUM62_MATRIX[i, j]
end

_aachar(c::Char) = uppercase(c)
_aachar(c::UInt8) = uppercase(Char(c))
_aachar(c::AbstractString) = (length(c) == 1 || throw(ArgumentError("expected a single amino acid")); uppercase(only(c)))
_aachar(c::Symbol) = _aachar(String(c))

"""
    RadicalCriteria(; grantham_min = 101, blosum62_max = nothing,
                    abs_delta_kd_min = nothing, combine = :any)

Definition of a *radical* amino acid substitution. Each criterion set to a
number is evaluated; `combine = :any` calls a substitution radical when at
least one criterion holds, `:all` when all hold.

The default (Grantham distance ≥ 101) follows the classes of Li, Wu & Luo
(1984): conservative 0–50, moderately conservative 51–100, moderately
radical 101–150, radical > 150. It is a convention, not a biological
constant.
"""
struct RadicalCriteria
    grantham_min::Union{Nothing,Float64}
    blosum62_max::Union{Nothing,Float64}
    abs_delta_kd_min::Union{Nothing,Float64}
    combine::Symbol
end

function RadicalCriteria(; grantham_min = 101, blosum62_max = nothing,
                         abs_delta_kd_min = nothing, combine::Symbol = :any)
    combine in (:any, :all) || throw(ArgumentError("combine must be :any or :all"))
    f(x) = x === nothing ? nothing : Float64(x)
    all(isnothing, (grantham_min, blosum62_max, abs_delta_kd_min)) &&
        throw(ArgumentError("RadicalCriteria needs at least one criterion"))
    RadicalCriteria(f(grantham_min), f(blosum62_max), f(abs_delta_kd_min), combine)
end

function Base.show(io::IO, r::RadicalCriteria)
    parts = String[]
    r.grantham_min === nothing || push!(parts, "Grantham ≥ $(r.grantham_min)")
    r.blosum62_max === nothing || push!(parts, "BLOSUM62 ≤ $(r.blosum62_max)")
    r.abs_delta_kd_min === nothing || push!(parts, "|ΔKD| ≥ $(r.abs_delta_kd_min)")
    print(io, join(parts, r.combine == :any ? " or " : " and "))
end

"""
    is_radical(a, b, criteria = RadicalCriteria()) -> Union{Bool,Missing}

Whether the substitution `a → b` between standard residues is radical.
"""
function is_radical(a, b, rc::RadicalCriteria = RadicalCriteria())
    g = grantham_distance(a, b)
    g === missing && return missing
    a == b && return false
    tests = Bool[]
    rc.grantham_min === nothing || push!(tests, g >= rc.grantham_min)
    rc.blosum62_max === nothing || push!(tests, blosum62_score(a, b) <= rc.blosum62_max)
    rc.abs_delta_kd_min === nothing ||
        push!(tests, abs(kd_value(_aachar(b)) - kd_value(_aachar(a))) >= rc.abs_delta_kd_min)
    rc.combine == :any ? any(tests) : all(tests)
end
