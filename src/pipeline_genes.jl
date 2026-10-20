# Per-gene preparation: gene extraction, de-duplication, reading frames,
# translation and per-sequence metrics. Every unique gene sequence is parsed
# and translated once; results are mapped back to the input rows.

mutable struct GeneData
    spec::GeneSpec
    row_unique::Vector{Int}              # unique gene sequence of each row (0 = none)
    row_status::Vector{Symbol}           # :ok, :truncated, :absent, :missing_sequence, :empty_sequence
    masks::Vector{Vector{UInt8}}         # unique gene segments (oriented)
    unique_rows::Vector{Vector{Int}}     # rows sharing each unique segment
    stopcounts::Vector{NTuple{3,Tuple{Int,Int}}}
    frame_status::Vector{Symbol}         # three-frame evaluation of each unique
    frame::Union{Int,Nothing}
    frame_source::Symbol                 # :user, :consensus, :unresolved
    frame_note::String
    records::Vector{Union{Missing,CodingRecord}}
    proteins::Vector{Union{Missing,String}}   # ungapped-by-column protein (gap codons kept as '-')
    hydro::Vector{Any}
end

genename(g::GeneData) = g.spec.name
nunique(g::GeneData) = length(g.masks)

"""
Deterministic threaded loop: `f(i)` for `i in 1:n`, each `i` exactly once.
Callers write to preallocated slots, so results do not depend on scheduling.
"""
function _tforeach(f, n::Int; threaded::Bool = true)
    if !threaded || Threads.nthreads() == 1 || n < 2
        for i in 1:n
            f(i)
        end
        return nothing
    end
    nchunks = min(n, 8 * Threads.nthreads())
    tasks = map(1:nchunks) do c
        Threads.@spawn for i in c:nchunks:n
            f(i)
        end
    end
    foreach(wait, tasks)
    nothing
end

function _parse_column(df::DataFrame, col)
    col === nothing && return nothing
    hasproperty(df, col) || throw(ArgumentError("column $(repr(col)) not found in the input DataFrame"))
    v = df[!, col]
    out = Vector{Union{Missing,ParsedSequence}}(undef, length(v))
    for i in eachindex(v)
        out[i] = parse_sequence(v[i])
    end
    out
end

function prepare_gene(spec::GeneSpec, df::DataFrame, main::Union{Nothing,Vector},
                      cfg::AssessmentConfig, colcache::Dict)
    n = nrow(df)
    src = if spec.column === nothing
        main === nothing && throw(ArgumentError(
            "gene $(spec.name): no sequence column; give `sequence_col` or a gene `column`"))
        main
    else
        get!(() -> _parse_column(df, spec.column), colcache, spec.column)
    end
    row_unique = zeros(Int, n)
    row_status = fill(:ok, n)
    index = Dict{Vector{UInt8},Int}()
    masks = Vector{Vector{UInt8}}()
    unique_rows = Vector{Vector{Int}}()
    for r in 1:n
        p = src[r]
        if p === missing
            row_status[r] = :missing_sequence
            continue
        elseif isempty(p.mask)
            row_status[r] = :empty_sequence
            continue
        end
        seg, st = extract_gene(p.mask, spec)
        row_status[r] = st
        st == :absent && continue
        if isempty(ungapped(seg))
            row_status[r] = :absent
            continue
        end
        u = get!(index, seg) do
            push!(masks, seg); push!(unique_rows, Int[])
            length(masks)
        end
        push!(unique_rows[u], r)
        row_unique[r] = u
    end
    code = genetic_code(cfg.genetic_code)
    nu = length(masks)
    stopcounts = Vector{NTuple{3,Tuple{Int,Int}}}(undef, nu)
    _tforeach(nu; threaded = cfg.threaded) do u
        stopcounts[u] = frame_stop_counts(masks[u], code; cfg.terminal_stop_allowed)
    end
    fstatus = [frame_status(sc) for sc in stopcounts]
    frame, source, note = _choose_frame(spec, stopcounts, fstatus, cfg)
    records = Vector{Union{Missing,CodingRecord}}(missing, nu)
    proteins = Vector{Union{Missing,String}}(missing, nu)
    hydro = Vector{Any}(undef, nu)
    _tforeach(nu; threaded = cfg.threaded) do u
        if frame === nothing
            hydro[u] = _missing_hydro("reading frame unresolved")
        else
            rec = read_frame(masks[u], frame, code; cfg.terminal_stop_allowed)
            records[u] = rec
            proteins[u] = rec.protein
            # Only standard residues: gaps, stops and X are skipped so that a
            # defect already reported by the coding test (an internal stop)
            # does not also break the windows and produce a second flag.
            hydro[u] = hydrophobicity_metrics(filter(is_standard_aa, rec.protein);
                                              window = cfg.hydrophobicity_window,
                                              window_threshold = cfg.window_threshold,
                                              hydrophobic_cutoff = cfg.hydrophobic_cutoff)
        end
    end
    GeneData(spec, row_unique, row_status, masks, unique_rows, stopcounts, fstatus, frame,
             source, note, records, proteins, hydro)
end

"""
Reading frame of a gene. A user frame is used as given. Otherwise the
frame is inferred only from a dataset consensus: among unique sequences
with exactly one open frame, one frame must be supported by at least
`frame_min_support` sequences and a share ≥ `frame_consensus_min`.
All sequences are then read in that frame (a sequence is never moved to
its own fewest-stop frame). Without consensus the frame stays unresolved.
"""
function _choose_frame(spec::GeneSpec, stopcounts, fstatus, cfg::AssessmentConfig)
    spec.frame !== nothing && return (spec.frame, :user, "frame given by the user")
    counts = zeros(Int, 3)
    for (sc, st) in zip(stopcounts, fstatus)
        st == :unique_open_frame || continue
        counts[findfirst(x -> x[1] == 0 && x[2] > 0, sc)] += 1
    end
    tot = sum(counts)
    nmulti = count(==(:multiple_open_frames), fstatus)
    nnone = count(==(:no_open_frame), fstatus)
    desc = "unique sequences with a single open frame: frame 1 = $(counts[1]), " *
           "frame 2 = $(counts[2]), frame 3 = $(counts[3]); several open frames: $nmulti; " *
           "no open frame: $nnone"
    if tot >= cfg.frame_min_support
        best = argmax(counts)
        share = counts[best] / tot
        if share >= cfg.frame_consensus_min
            return (best, :consensus,
                    "frame inferred from dataset consensus ($(round(share; digits = 3)) of $tot); " * desc)
        end
    end
    (nothing, :unresolved, "no frame consensus; " * desc)
end

_open_frames(sc) = join([string(f) for f in 1:3 if sc[f][1] == 0 && sc[f][2] > 0], ",")

function validation_table(g::GeneData, ids, invalid_counts)
    n = length(g.row_unique)
    T = Union{Missing,Int}; F = Union{Missing,Float64}
    cols = (frame = T[], length_nt = T[], n_columns = T[], n_ambiguous_nt = T[], gc = F[],
            gc1 = F[], gc2 = F[], gc3 = F[], leading_bases = T[], n_codons = T[],
            trailing_bases = T[], n_internal_stops = T[], terminal_stop = Union{Missing,Bool}[],
            n_ambiguous_codons = T[], n_unresolved_codons = T[], n_possible_stops = T[],
            n_partial_gap_codons = T[], n_gap_codons = T[])
    stopcodons = Union{Missing,String}[]; openf = Union{Missing,String}[]
    fstat = Union{Missing,Symbol}[]; prot = Union{Missing,String}[]
    for r in 1:n
        u = g.row_unique[r]
        rec = u == 0 ? missing : g.records[u]
        if rec === missing
            for c in values(cols)
                push!(c, missing)
            end
            if u != 0       # sequence present, frame unresolved: frame-free fields
                m = g.masks[u]
                cols.length_nt[end] = count(!=(GAP_MASK), m)
                cols.n_columns[end] = length(m)
                cols.n_ambiguous_nt[end] = count(is_ambiguous, m)
                cols.gc[end] = _gc(m)
            end
            push!(stopcodons, missing); push!(prot, missing)
        else
            for f in keys(cols)
                f == :frame && (push!(cols.frame, rec.frame); continue)
                push!(cols[f], getfield(rec, f))
            end
            push!(stopcodons, join(rec.internal_stop_positions, ","))
            push!(prot, rec.protein)
        end
        push!(openf, u == 0 ? missing : _open_frames(g.stopcounts[u]))
        push!(fstat, u == 0 ? missing : g.frame_status[u])
    end
    df = DataFrame(row = 1:n, id = ids, gene = fill(genename(g), n), status = copy(g.row_status),
                   frame_source = fill(g.frame_source, n))
    for (k, v) in pairs(cols)
        df[!, k] = v
    end
    df.internal_stop_codons = stopcodons
    df.open_frames = openf
    df.three_frame_status = fstat
    df.n_invalid_symbols = invalid_counts
    df.protein = prot
    df
end

function hydrophobicity_rows(g::GeneData, ids)
    n = length(g.row_unique)
    rows = [g.row_unique[r] == 0 ? _missing_hydro(string(g.row_status[r])) : g.hydro[g.row_unique[r]]
            for r in 1:n]
    df = DataFrame(rows)
    insertcols!(df, 1, :row => 1:n, :id => ids, :gene => fill(genename(g), n))
    df
end

function composition_rows(g::GeneData, ids)
    n = length(g.row_unique)
    comp = Dict(u => (g.proteins[u] === missing ? fill(missing, 20) : aa_composition(g.proteins[u]))
                for u in 1:nunique(g))
    df = DataFrame(row = 1:n, id = ids, gene = fill(genename(g), n))
    for (k, aa) in enumerate(AMINO_ACIDS)
        df[!, Symbol(aa)] = Union{Missing,Float64}[g.row_unique[r] == 0 ? missing :
                                                   comp[g.row_unique[r]][k] for r in 1:n]
    end
    df
end

"""Analysis state shared by the comparative steps."""
struct AnalysisContext
    input::DataFrame
    ids::Vector
    cfg::AssessmentConfig
    genes::Vector{GeneData}
    main::Union{Nothing,Vector}
    abundance::Union{Nothing,Vector{Union{Missing,Float64}}}
    colcache::Dict{Symbol,Any}
end
