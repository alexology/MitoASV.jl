# Sequence parsing, gene specifications and gene extraction.

"""
    ParsedSequence

Analysis copy of one input sequence. `mask` has one IUPAC mask per input
character. Characters that are not IUPAC nucleotide codes or gaps
(`-`, `.`) are counted in `n_invalid` and treated as `N` (unknown) during
analysis. The input itself is never changed.
"""
struct ParsedSequence
    mask::Vector{UInt8}
    n_invalid::Int
    n_rna::Int          # number of `U`/`u` symbols (read as T)
end

_as_string(x::AbstractString) = x
_as_string(x::AbstractVector{<:AbstractChar}) = String(x)
_as_string(x) = string(x)       # BioSequences.LongDNA/LongRNA, LongSequence, etc.

"""
    parse_sequence(x) -> ParsedSequence or missing

Parse a nucleotide sequence given as a string, a character vector or any
object whose `string` representation is a nucleotide string (for example
`BioSequences.LongDNA`). `missing` and `nothing` give `missing`.
"""
parse_sequence(::Missing) = missing
parse_sequence(::Nothing) = missing
function parse_sequence(x)
    s = _as_string(x)
    mask = Vector{UInt8}(undef, length(s))
    ninv = 0
    nrna = 0
    i = 0
    for c in s
        i += 1
        m = isascii(c) ? NT_MASK[UInt8(c) + 1] : INVALID_MASK
        if m == INVALID_MASK
            ninv += 1
            m = 0x0f
        end
        nrna += (c == 'U' || c == 'u')
        mask[i] = m
    end
    ParsedSequence(mask, ninv, nrna)
end

# --- Gene specifications ----------------------------------------------------

"""
    GeneSpec(name; start = 1, stop = nothing, frame = nothing, strand = :+, column = nothing)

Location of one protein-coding gene.

- `start`, `stop`: **1-based inclusive** positions on the input sequence as
  given (forward orientation; alignment columns if the input is aligned).
  `stop = nothing` means "to the end of the sequence".
- `frame`: 1, 2 or 3, the position of the first base of the first complete
  codon, counted from the 5' end of the *oriented* gene segment. `nothing`
  means unknown; the frame is then inferred from a dataset consensus (see
  `assess_asvs`) or left unresolved.
- `strand`: `:+` (as given) or `:-` (the gene is the reverse complement of
  the `start:stop` segment).
- `column`: name of a column that already holds the extracted gene
  sequence. When set, `start`/`stop` refer to that column's sequences.

A `NamedTuple` with these fields is accepted wherever a `GeneSpec` is
expected, e.g. `(start = 1, stop = 900, frame = 1)`.
"""
struct GeneSpec
    name::Symbol
    start::Int
    stop::Union{Int,Nothing}
    frame::Union{Int,Nothing}
    strand::Symbol
    column::Union{Symbol,Nothing}
end

function GeneSpec(name; start::Integer = 1, stop = nothing, frame = nothing, strand = :+,
                  column = nothing)
    start >= 1 || throw(ArgumentError("gene $name: start must be ≥ 1 (1-based coordinates)"))
    stop === nothing || stop >= start ||
        throw(ArgumentError("gene $name: stop ($stop) must be ≥ start ($start)"))
    frame === nothing || frame in 1:3 ||
        throw(ArgumentError("gene $name: frame must be 1, 2, 3 or nothing"))
    st = Symbol(strand)
    st = st in (:+, :forward, :plus) ? :+ : st in (:-, :reverse, :minus) ? :- :
         throw(ArgumentError("gene $name: strand must be :+ or :-"))
    GeneSpec(Symbol(name), Int(start), stop === nothing ? nothing : Int(stop),
             frame === nothing ? nothing : Int(frame), st,
             column === nothing ? nothing : Symbol(column))
end
GeneSpec(name, nt::NamedTuple) = GeneSpec(name; nt...)
GeneSpec(name, column::Union{Symbol,AbstractString}) = GeneSpec(name; column = Symbol(column))
GeneSpec(name, g::GeneSpec) = GeneSpec(Symbol(name), g.start, g.stop, g.frame, g.strand, g.column)

"""
Extract the gene segment from a parsed sequence. Returns
`(segment, status)` with status `:ok`, `:truncated` (the sequence ends
before `stop`) or `:absent` (the sequence ends before `start`).
"""
function extract_gene(mask::Vector{UInt8}, spec::GeneSpec)
    n = length(mask)
    spec.start > n && return (UInt8[], :absent)
    stop = spec.stop === nothing ? n : spec.stop
    status = stop > n ? :truncated : :ok
    seg = mask[spec.start:min(stop, n)]
    if spec.strand == :-
        reverse!(seg)
        seg .= complement_mask.(seg)
    end
    seg, status
end

ungapped(mask::AbstractVector{UInt8}) = filter(!=(GAP_MASK), mask)

"""
    reverse_complement(seq::AbstractString) -> String

Reverse complement of an IUPAC nucleotide string (gaps are kept).
"""
function reverse_complement(seq::AbstractString)
    p = parse_sequence(seq)
    mask_string(complement_mask.(reverse(p.mask)))
end

"""
    gc_content(seq) -> Union{Float64,Missing}

Fraction of G+C among bases whose GC status is known (A, C, G, T, S, W).
Other ambiguity codes and gaps are excluded; `missing` if none remain.
"""
gc_content(seq) = (p = parse_sequence(seq); p === missing ? missing : _gc(p.mask))

function _gc(mask::AbstractVector{UInt8}, idx = eachindex(mask))
    gc = 0
    tot = 0
    for i in idx
        m = mask[i]
        if is_gc(m)
            gc += 1; tot += 1
        elseif is_at(m)
            tot += 1
        end
    end
    tot == 0 ? missing : gc / tot
end
