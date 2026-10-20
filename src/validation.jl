# Coding sequence validation: reading frames, translation, stop codons,
# ambiguous codons, base composition.

"""
    CodingRecord

Result of reading one (gene) sequence in one frame. Codon positions refer
to the oriented gene segment. Gap columns are kept so that codon `k` is the
same alignment column for every sequence of an aligned input.
"""
struct CodingRecord
    frame::Int
    length_nt::Int               # non-gap positions in the segment
    n_columns::Int               # positions including gaps
    n_ambiguous_nt::Int
    gc::Union{Missing,Float64}
    gc1::Union{Missing,Float64}
    gc2::Union{Missing,Float64}
    gc3::Union{Missing,Float64}
    leading_bases::Int           # columns before the first complete codon
    n_codons::Int                # complete codons (including full gap codons)
    trailing_bases::Int          # columns after the last complete codon (0, 1, 2)
    protein::String              # one character per complete codon ('-' = gap codon)
    n_internal_stops::Int
    internal_stop_positions::Vector{Int}   # codon indices
    terminal_stop::Bool
    n_ambiguous_codons::Int      # codons with ≥1 ambiguous base
    n_unresolved_codons::Int     # ambiguous codons translated as 'X'
    n_possible_stops::Int        # ambiguous codons that may be stops
    n_partial_gap_codons::Int    # codons with 1 or 2 gap columns
    n_gap_codons::Int
end

"""
Read `mask` in `frame` with `code`. With `terminal_stop_allowed = true` a
stop in the last non-gap codon, with no non-gap bases after it, is reported
as `terminal_stop` and not as internal.
"""
function read_frame(mask::AbstractVector{UInt8}, frame::Int, code::GeneticCode;
                    terminal_stop_allowed::Bool = true)
    n = length(mask)
    lead = min(frame - 1, n)
    ncod = max(n - lead, 0) ÷ 3
    trail = n - lead - 3 * ncod
    prot = Vector{UInt8}(undef, ncod)
    namb_cod = 0; nunres = 0; nmaybe = 0; npgap = 0; ngap = 0
    stops = Int[]
    @inbounds for k in 1:ncod
        p = lead + 3 * (k - 1)
        m1, m2, m3 = mask[p + 1], mask[p + 2], mask[p + 3]
        aa = translate_codon(code, m1, m2, m3)
        prot[k] = aa
        if aa == UInt8('-')
            ngap += 1
        elseif m1 == GAP_MASK || m2 == GAP_MASK || m3 == GAP_MASK
            npgap += 1
        elseif is_ambiguous(m1) || is_ambiguous(m2) || is_ambiguous(m3)
            namb_cod += 1
            aa == UInt8('X') && (nunres += 1)
            code.maybe_stop[_lut_index(m1, m2, m3)] && (nmaybe += 1)
        end
        aa == UInt8('*') && push!(stops, k)
    end
    # Terminal stop: the last non-gap codon, with no non-gap bases after it.
    terminal = false
    if terminal_stop_allowed && !isempty(stops)
        last_nongap = findlast(!=(UInt8('-')), prot)
        tail_bases = count(i -> mask[i] != GAP_MASK, (lead + 3 * ncod + 1):n)
        if last_nongap !== nothing && stops[end] == last_nongap && tail_bases == 0
            terminal = true
            pop!(stops)
        end
    end
    nnt = count(!=(GAP_MASK), mask)
    namb = count(is_ambiguous, mask)
    pos(r) = (lead + r):3:(lead + 3 * ncod)
    CodingRecord(frame, nnt, n, namb, _gc(mask), _gc(mask, pos(1)), _gc(mask, pos(2)),
                 _gc(mask, pos(3)), lead, ncod, trail, String(prot), length(stops), stops,
                 terminal, namb_cod, nunres, nmaybe, npgap, ngap)
end

"""
Count internal stops in each of the three forward frames.
Returns a 3-tuple of `(n_internal_stops, n_codons)`.
"""
function frame_stop_counts(mask::AbstractVector{UInt8}, code::GeneticCode;
                           terminal_stop_allowed::Bool = true)
    ntuple(f -> begin
        r = read_frame(mask, f, code; terminal_stop_allowed)
        (r.n_internal_stops, r.n_codons)
    end, 3)
end

"""
    frame_candidates(seq; genetic_code = 5, strand = :+, terminal_stop_allowed = true) -> DataFrame

Evaluate the three reading frames of one sequence: number of complete
codons, internal stops and whether the frame is open (no internal stop).
No frame is selected: several open frames, or none, are reported as such.
"""
function frame_candidates(seq; genetic_code = 5, strand = :+, terminal_stop_allowed = true)
    p = parse_sequence(seq)
    p === missing && throw(ArgumentError("sequence is missing"))
    mask = Symbol(strand) in (:-, :reverse, :minus) ? complement_mask.(reverse(p.mask)) : p.mask
    code = MitoASV.genetic_code(genetic_code)
    rows = [read_frame(mask, f, code; terminal_stop_allowed) for f in 1:3]
    DataFrame(frame = 1:3, n_codons = [r.n_codons for r in rows],
              n_internal_stops = [r.n_internal_stops for r in rows],
              open = [r.n_internal_stops == 0 && r.n_codons > 0 for r in rows],
              trailing_bases = [r.trailing_bases for r in rows])
end

"""
    frame_status(stopcounts) -> Symbol

`:unique_open_frame`, `:multiple_open_frames` or `:no_open_frame`.
"""
function frame_status(sc)
    nopen = count(x -> x[1] == 0 && x[2] > 0, sc)
    nopen == 1 ? :unique_open_frame : nopen > 1 ? :multiple_open_frames : :no_open_frame
end

"""
    translate_sequence(seq; genetic_code = 5, frame = 1, strand = :+) -> String

Translate a nucleotide sequence from `frame` (1, 2 or 3). Bases before the
frame and an incomplete terminal codon are not translated. Ambiguous codons
become the amino acid shared by all possible codons, or `X`. Full gap
codons become `-`, partial gap codons `X`. Stops are `*`.
"""
function translate_sequence(seq; genetic_code = 5, frame::Integer = 1, strand = :+)
    frame in 1:3 || throw(ArgumentError("frame must be 1, 2 or 3"))
    p = parse_sequence(seq)
    p === missing && return missing
    mask = Symbol(strand) in (:-, :reverse, :minus) ? complement_mask.(reverse(p.mask)) : p.mask
    read_frame(mask, Int(frame), MitoASV.genetic_code(genetic_code)).protein
end

"""
    coding_validation(seq; frame = nothing, genetic_code = 5, strand = :+,
                      terminal_stop_allowed = true) -> NamedTuple

Validate one coding sequence. With `frame = nothing` the three frames are
evaluated and reported (`frame_table`, `frame_status`) but **no frame is
chosen** and the frame-specific fields are `missing`; give `frame` to get
them. An incomplete terminal codon is reported in `trailing_bases` and is
not treated as a frameshift.
"""
function coding_validation(seq; frame = nothing, genetic_code = 5, strand = :+,
                           terminal_stop_allowed = true)
    p = parse_sequence(seq)
    p === missing && throw(ArgumentError("sequence is missing"))
    code = MitoASV.genetic_code(genetic_code)
    mask = Symbol(strand) in (:-, :reverse, :minus) ? complement_mask.(reverse(p.mask)) : p.mask
    sc = frame_stop_counts(mask, code; terminal_stop_allowed)
    status = frame_status(sc)
    ft = DataFrame(frame = 1:3, n_codons = [s[2] for s in sc],
                   n_internal_stops = [s[1] for s in sc],
                   open = [s[1] == 0 && s[2] > 0 for s in sc])
    base = (length_nt = count(!=(GAP_MASK), mask), n_invalid_symbols = p.n_invalid,
            n_ambiguous_nt = count(is_ambiguous, mask), gc = _gc(mask),
            frame_status = status, frame_table = ft)
    if frame === nothing
        return merge(base, (frame = missing, record = missing))
    end
    frame in 1:3 || throw(ArgumentError("frame must be 1, 2 or 3"))
    r = read_frame(mask, Int(frame), code; terminal_stop_allowed)
    merge(base, (frame = Int(frame), record = r, protein = r.protein,
                 n_internal_stops = r.n_internal_stops, terminal_stop = r.terminal_stop,
                 trailing_bases = r.trailing_bases, gc1 = r.gc1, gc2 = r.gc2, gc3 = r.gc3))
end

const AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

"""
    aa_composition(protein) -> Vector{Float64}

Fractions of the 20 standard amino acids (order `$(AMINO_ACIDS)`) among
standard residues; `X`, `*`, `-` and other symbols are excluded. All
`missing` if the protein has no standard residue.
"""
function aa_composition(protein::AbstractString)
    counts = zeros(Int, 20)
    for c in protein
        k = findfirst(==(c), AMINO_ACIDS)
        k === nothing || (counts[k] += 1)
    end
    tot = sum(counts)
    out = Vector{Union{Missing,Float64}}(undef, 20)
    tot == 0 ? fill!(out, missing) : (out .= counts ./ tot)
end
