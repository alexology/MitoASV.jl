# Reference-free conservation estimated from the ASVs themselves.
#
# Conservation is estimated within comparison groups of aligned
# (equal-length) proteins. Each sequence is scored against frequencies
# computed WITHOUT it (leave-one-out) to avoid self-influence. ASVs are not
# independent evolutionary observations: closely related or duplicated
# sequences inflate apparent conservation. Identical sequences are collapsed
# before this step in `assess_asvs`; optional Henikoff weights further
# down-weight redundant sequences.

const LOG2_20 = log2(20.0)

struct ConservationCore
    L::Int
    W::Matrix{Float64}      # L × 20 weighted residue counts
    R::Matrix{Int}          # L × 20 raw counts
    wsum::Vector{Float64}
    rsum::Vector{Int}
end

function conservation_core(prots::AbstractVector, w::AbstractVector{Float64})
    L = length(prots[1])
    W = zeros(L, 20); R = zeros(Int, L, 20)
    for (s, wi) in zip(prots, w)
        @inbounds for p in 1:L
            k = aa_index(s[p])
            k == 0 && continue
            W[p, k] += wi; R[p, k] += 1
        end
    end
    ConservationCore(L, W, R, vec(sum(W; dims = 2)), vec(sum(R; dims = 2)))
end

"""
    henikoff_weights(proteins) -> Vector{Float64}

Position-based sequence weights (Henikoff & Henikoff 1994), scaled to mean
1. Redundant sequences share weight; a sequence with no standard residue
gets weight 1.
"""
function henikoff_weights(prots::AbstractVector)
    n = length(prots); L = length(prots[1])
    w = zeros(n); npos = 0
    cnt = zeros(Int, 20)
    for p in 1:L
        fill!(cnt, 0)
        for s in prots
            k = aa_index(s[p]); k == 0 || (cnt[k] += 1)
        end
        r = count(>(0), cnt)
        r == 0 && continue
        npos += 1
        for (i, s) in enumerate(prots)
            k = aa_index(s[p]); k == 0 || (w[i] += 1 / (r * cnt[k]))
        end
    end
    npos == 0 && return ones(n)
    for i in 1:n
        w[i] == 0 && (w[i] = NaN)
    end
    m = mean(filter(!isnan, w))
    [isnan(x) ? 1.0 : x / m for x in w]
end

n_effective(w) = (s = sum(w); s2 = sum(abs2, w); s2 == 0 ? 0.0 : s^2 / s2)

"""
Leave-one-out frequencies at position `p` without residue index `k0`
(weight `w0`). Returns `(consensus, consensus_freq, tie, conservation,
freq_of(k), informative)`.
"""
@inline function _loo_position(core::ConservationCore, p::Int, k0::Int, w0::Float64, min_inf::Int)
    tot = core.wsum[p] - (k0 == 0 ? 0.0 : w0)
    raw = core.rsum[p] - (k0 == 0 ? 0 : 1)
    (raw < min_inf || tot <= 0) && return (0, 0.0, false, NaN, tot, false)
    best = 0; bestv = -1.0; tie = false; H = 0.0
    @inbounds for k in 1:20
        v = core.W[p, k] - (k == k0 ? w0 : 0.0)
        v < 1e-12 && continue
        f = v / tot
        H -= f * log2(f)
        if v > bestv + 1e-9
            best = k; bestv = v; tie = false
        elseif abs(v - bestv) <= 1e-9
            tie = true
        end
    end
    (best, bestv / tot, tie, 1 - H / LOG2_20, tot, true)
end

"""
Leave-one-out amino acid and conservation scores of one protein `s`
(bytes) with weight `w0` against `core`.
"""
function loo_protein_scores(core::ConservationCore, s, w0::Float64, min_inf::Int,
                            rare_threshold::Float64, radical::RadicalCriteria)
    ncomp = 0; nsub = 0; gsum = 0; gmax = 0; bsum = 0; bmin = typemax(Int)
    dkd = 0.0; adkd = 0.0; nrad = 0; nstop = 0; ngap = 0; nnocons = 0
    nscored = 0; csum = 0.0; nunseen = 0; nrare = 0
    cwg = 0.0; cwgmax = 0.0; cwkd = 0.0; nradinv = 0
    @inbounds for p in 1:core.L
        c = s[p]
        k0 = aa_index(c)
        cons, cf, tie, cons_c, tot, inf = _loo_position(core, p, k0, w0, min_inf)
        if !inf
            nnocons += 1
            continue
        end
        nscored += 1; csum += cons_c
        if k0 != 0
            f = (core.W[p, k0] - w0) / tot
            f < 1e-12 && (nunseen += 1)
            f < rare_threshold && (nrare += 1)
        end
        if tie
            nnocons += 1
            continue
        end
        if c == UInt8('*')
            nstop += 1; continue
        elseif c == UInt8('-')
            ngap += 1; continue
        elseif k0 == 0
            continue
        end
        ncomp += 1
        k0 == cons && continue
        nsub += 1
        a = AMINO_ACIDS[cons]; b = Char(c)
        g = GRANTHAM_MATRIX[cons, k0]; bl = blosum62_score(a, b)
        d = KD_VEC[k0] - KD_VEC[cons]
        gsum += g; gmax = max(gmax, g); bsum += bl; bmin = min(bmin, bl)
        dkd += d; adkd += abs(d)
        rad = is_radical(a, b, radical)
        rad && (nrad += 1)
        cwg += cons_c * g; cwgmax = max(cwgmax, cons_c * g); cwkd += cons_c * abs(d)
        rad && cf >= 1 - 1e-12 && (nradinv += 1)
    end
    aa = (n_compared = ncomp, n_substitutions = nsub,
          prop_different = ncomp == 0 ? missing : nsub / ncomp,
          grantham_sum = ncomp == 0 ? missing : gsum,
          grantham_mean = nsub == 0 ? missing : gsum / nsub,
          grantham_max = nsub == 0 ? missing : gmax,
          blosum62_mean = nsub == 0 ? missing : bsum / nsub,
          blosum62_min = nsub == 0 ? missing : bmin,
          delta_kd_sum = ncomp == 0 ? missing : dkd,
          abs_delta_kd_sum = ncomp == 0 ? missing : adkd,
          n_radical = ncomp == 0 ? missing : nrad,
          n_stop_vs_consensus = nstop, n_gap_vs_consensus = ngap,
          n_positions_without_consensus = nnocons)
    cons = (n_scored_positions = nscored,
            mean_conservation = nscored == 0 ? missing : csum / nscored,
            n_unseen_residues = nscored == 0 ? missing : nunseen,
            n_rare_residues = nscored == 0 ? missing : nrare,
            cw_grantham_sum = ncomp == 0 ? missing : cwg,
            cw_grantham_max = ncomp == 0 ? missing : cwgmax,
            cw_abs_delta_kd_sum = ncomp == 0 ? missing : cwkd,
            n_radical_invariant = ncomp == 0 ? missing : nradinv)
    aa, cons
end

"""
    conservation_analysis(sequences; weights = nothing, min_sequences = 5,
                          rare_threshold = 0.05, radical = RadicalCriteria()) -> NamedTuple

Conservation of aligned (equal-length) protein sequences estimated from the
sequences themselves; `missing` entries are skipped.

- `weights = nothing`: unweighted (each sequence counts once). Other
  options: `:henikoff` (position-based weights that down-weight redundant
  sequences), `:unique` (identical sequences collapsed to one), or a
  numeric vector (e.g. abundances; reads are not independent
  evolutionary observations, so interpret with care).
- Positions need at least `min_sequences - 1` informative residues in the
  leave-one-out set, and the analysis needs at least `min_sequences`
  sequences; otherwise values are `missing` and `reason` explains why.

Returns
- `positions`: per position `n_informative`, `consensus` (`missing` on
  ties), `consensus_freq`, `entropy` (bits, 20 residues), `conservation`
  (1 − entropy / log2 20) and `freq_A` … `freq_Y`, using all sequences;
- `sequences`: per sequence leave-one-out scores: substitutions against the
  leave-one-out consensus (Grantham, BLOSUM62, hydropathy change, radical
  changes), unseen and rare residues, conservation-weighted Grantham and
  hydropathy changes, and radical substitutions at sites invariant among
  the other sequences;
- `n_sequences`, `n_effective` (Kish effective number from the weights),
  `weights`, `reason`.
"""
function conservation_analysis(sequences::AbstractVector; weights = nothing, min_sequences::Integer = 5,
                               rare_threshold::Real = 0.05,
                               radical::RadicalCriteria = RadicalCriteria())
    keep = findall(!ismissing, sequences)
    prots = [Vector{UInt8}(codeunits(uppercase(String(sequences[i])))) for i in keep]
    n = length(prots)
    n > 0 && length(unique(length.(prots))) == 1 || n == 0 || throw(ArgumentError(
        "sequences differ in length; conservation needs aligned (equal-length) proteins"))
    w = if weights === nothing
        ones(n)
    elseif weights === :henikoff
        n == 0 ? Float64[] : henikoff_weights(prots)
    elseif weights === :unique
        ones(n)              # replaced below
    elseif weights isa AbstractVector
        length(weights) == length(sequences) ||
            throw(ArgumentError("weights must have one value per sequence"))
        ww = Float64[Float64(weights[i]) for i in keep]
        all(x -> x >= 0 && isfinite(x), ww) || throw(ArgumentError("weights must be finite and ≥ 0"))
        ww
    else
        throw(ArgumentError("weights must be nothing, :henikoff, :unique or a vector"))
    end
    if weights === :unique   # first occurrence of each distinct sequence gets weight 1
        seen = Set{Vector{UInt8}}()
        w = [p in seen ? 0.0 : (push!(seen, p); 1.0) for p in prots]
    end
    neff = n == 0 ? 0.0 : n_effective(w)
    if n < min_sequences
        reason = "fewer than $min_sequences sequences ($n)"
        return (positions = DataFrame(), sequences = _empty_seq_scores(sequences, reason),
                n_sequences = n, n_effective = neff, weights = w, reason = reason)
    end
    core = conservation_core(prots, w)
    min_inf = min_sequences - 1
    # positions (all sequences)
    pos = DataFrame(position = 1:core.L, n_informative = core.rsum)
    consv = Union{Missing,Char}[]; cf = Union{Missing,Float64}[]; ent = Union{Missing,Float64}[]
    conserv = Union{Missing,Float64}[]; ties = Bool[]
    for p in 1:core.L
        c, f, tie, cc, _, inf = _loo_position(core, p, 0, 0.0, min_sequences)
        if !inf
            push!(consv, missing); push!(cf, missing); push!(ent, missing); push!(conserv, missing)
            push!(ties, false)
        else
            push!(consv, tie ? missing : AMINO_ACIDS[c]); push!(cf, f)
            push!(conserv, cc); push!(ent, (1 - cc) * LOG2_20); push!(ties, tie)
        end
    end
    pos.consensus = consv; pos.consensus_freq = cf; pos.tie = ties
    pos.entropy = ent; pos.conservation = conserv
    for (k, a) in enumerate(AMINO_ACIDS)
        pos[!, Symbol("freq_", a)] = [core.wsum[p] > 0 ? core.W[p, k] / core.wsum[p] : missing
                                      for p in 1:core.L]
    end
    rows = map(eachindex(prots)) do i
        aa, cs = loo_protein_scores(core, prots[i], w[i], min_inf, Float64(rare_threshold), radical)
        merge(aa, cs)
    end
    seqdf = DataFrame(rows)
    insertcols!(seqdf, 1, :index => keep, :weight => w)
    seqdf.reason = fill("", n)
    out = _empty_seq_scores(sequences, "missing sequence")
    (positions = pos, sequences = _merge_scores(out, seqdf), n_sequences = n,
     n_effective = neff, weights = w, reason = "")
end

function _empty_seq_scores(sequences, reason)
    DataFrame(index = collect(eachindex(sequences)), reason = fill(reason, length(sequences)))
end

function _merge_scores(base, seqdf)
    out = leftjoin(select(base, :index), seqdf; on = :index)
    sort!(out, :index)
    out.reason = coalesce.(out.reason, "missing sequence")
    out
end
