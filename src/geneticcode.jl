# Nucleotide encoding and NCBI genetic codes.
#
# Sequences are stored internally as vectors of IUPAC bit masks
# (A = 1, C = 2, G = 4, T/U = 8; ambiguity codes are unions; gap = 16).
# The original input is never modified: the mask vector is an analysis copy.

const GAP_MASK = 0x10
const INVALID_MASK = 0xff

const NT_MASK = let t = fill(INVALID_MASK, 256)
    for (c, m) in (('A', 1), ('C', 2), ('G', 4), ('T', 8), ('U', 8), ('R', 5),
                   ('Y', 10), ('S', 6), ('W', 9), ('K', 12), ('M', 3), ('B', 14),
                   ('D', 13), ('H', 11), ('V', 7), ('N', 15))
        t[UInt8(c) + 1] = m
        t[UInt8(lowercase(c)) + 1] = m
    end
    t[UInt8('-') + 1] = GAP_MASK
    t[UInt8('.') + 1] = GAP_MASK
    t
end

# Display character for each mask value 0..16.
const MASK_CHAR = ('?', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K',
                   'D', 'B', 'N', '-')

@inline is_base(m::UInt8) = m != GAP_MASK && m != 0x00
@inline is_unambiguous(m::UInt8) = m == 0x01 || m == 0x02 || m == 0x04 || m == 0x08
@inline is_ambiguous(m::UInt8) = is_base(m) && !is_unambiguous(m)
@inline is_gc(m::UInt8) = m == 0x02 || m == 0x04 || m == 0x06          # C, G, S
@inline is_at(m::UInt8) = m == 0x01 || m == 0x08 || m == 0x09          # A, T, W

@inline function complement_mask(m::UInt8)
    m == GAP_MASK && return m
    ((m & 0x01) << 3) | ((m & 0x08) >> 3) | ((m & 0x02) << 1) | ((m & 0x04) >> 1)
end

mask_string(mask::AbstractVector{UInt8}) = String([MASK_CHAR[Int(m) + 1] for m in mask])

# --- Genetic codes ----------------------------------------------------------

# NCBI translation tables in TCAG order (first base varies slowest).
# Source: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
const NCBI_TABLES = Dict{Int,Tuple{String,String}}(
    1 => ("Standard",
          "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    2 => ("Vertebrate Mitochondrial",
          "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"),
    4 => ("Mold, Protozoan, and Coelenterate Mitochondrial; Mycoplasma; Spiroplasma",
          "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"),
    5 => ("Invertebrate Mitochondrial",
          "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG"),
    9 => ("Echinoderm and Flatworm Mitochondrial",
          "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG"),
)

"""
    GeneticCode

An NCBI translation table. Construct with [`genetic_code`](@ref).
`table` holds the 64 amino acids in NCBI TCAG order. `lut` translates any
combination of three IUPAC masks: an ambiguous codon gets an amino acid
only when every possible codon gives the same one, otherwise `'X'`.
"""
struct GeneticCode
    id::Int
    name::String
    table::String
    lut::Vector{UInt8}          # 4096 entries indexed by (m1, m2, m3) masks 0..15
    maybe_stop::BitVector       # ambiguous codon that could be a stop but is not always one
end

Base.show(io::IO, c::GeneticCode) = print(io, "GeneticCode(", c.id, ", \"", c.name, "\")")

@inline _tcag(b::Integer) = b == 8 ? 0 : b == 2 ? 1 : b == 1 ? 2 : 3
_bits(m::Integer) = [b for b in (1, 2, 4, 8) if m & b != 0]
@inline _lut_index(m1::UInt8, m2::UInt8, m3::UInt8) =
    (Int(m1) << 8) | (Int(m2) << 4) | Int(m3) + 1

function _build_code(id::Int)
    name, table = NCBI_TABLES[id]
    lut = fill(UInt8('X'), 4096)
    maybe = falses(4096)
    for m1 in 1:15, m2 in 1:15, m3 in 1:15
        aas = Set{UInt8}()
        for b1 in _bits(m1), b2 in _bits(m2), b3 in _bits(m3)
            push!(aas, UInt8(table[16 * _tcag(b1) + 4 * _tcag(b2) + _tcag(b3) + 1]))
        end
        k = _lut_index(UInt8(m1), UInt8(m2), UInt8(m3))
        lut[k] = length(aas) == 1 ? first(aas) : UInt8('X')
        maybe[k] = length(aas) > 1 && UInt8('*') in aas
    end
    GeneticCode(id, name, table, lut, maybe)
end

const GENETIC_CODES = Dict{Int,GeneticCode}(id => _build_code(id) for id in keys(NCBI_TABLES))

"""
    genetic_code(id::Integer = 5) -> GeneticCode

Return NCBI translation table `id`. Supported: 1 (standard), 2 (vertebrate
mitochondrial), 4, 5 (invertebrate mitochondrial, the default for insects)
and 9.
"""
function genetic_code(id::Integer = 5)
    haskey(GENETIC_CODES, id) || throw(ArgumentError(
        "unsupported genetic code $id; supported NCBI tables: $(sort(collect(keys(GENETIC_CODES))))"))
    GENETIC_CODES[id]
end
genetic_code(c::GeneticCode) = c

"""
Translate one codon given as three masks. Returns the amino acid byte:
`'-'` for a full gap codon, `'X'` for a partial gap codon or an
unresolvable ambiguous codon.
"""
@inline function translate_codon(code::GeneticCode, m1::UInt8, m2::UInt8, m3::UInt8)
    if m1 == GAP_MASK || m2 == GAP_MASK || m3 == GAP_MASK
        return (m1 == GAP_MASK && m2 == GAP_MASK && m3 == GAP_MASK) ? UInt8('-') : UInt8('X')
    end
    code.lut[_lut_index(m1, m2, m3)]
end

"""
    translate_codon(codon::AbstractString; genetic_code = 5) -> Char

Translate a single codon such as `"TGA"` (tryptophan in table 5).
"""
function translate_codon(codon::AbstractString; genetic_code = 5)
    length(codon) == 3 || throw(ArgumentError("a codon must have exactly 3 nucleotides"))
    m = parse_sequence(codon).mask
    any(==(INVALID_MASK), m) && return 'X'
    Char(translate_codon(MitoASV.genetic_code(genetic_code), m[1], m[2], m[3]))
end

# Unambiguous codon index 1..64 in TCAG order, or 0 if ambiguous/gap.
@inline function codon_index(m1::UInt8, m2::UInt8, m3::UInt8)
    (is_unambiguous(m1) && is_unambiguous(m2) && is_unambiguous(m3)) || return 0
    16 * _tcag(m1) + 4 * _tcag(m2) + _tcag(m3) + 1
end

const TCAG = ('T', 'C', 'A', 'G')
codon_string(idx::Integer) = string(TCAG[(idx - 1) ÷ 16 + 1], TCAG[((idx - 1) ÷ 4) % 4 + 1],
                                    TCAG[(idx - 1) % 4 + 1])
