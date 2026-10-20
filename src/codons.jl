# Codon-level comparisons (Nei & Gojobori 1986), codon usage and RSCU.

"""
Per-code tables for unambiguous codons 1..64 (TCAG order):
`aa[c]`, synonymous sites `ssites[c]`, and the Nei-Gojobori pathway
average of synonymous / nonsynonymous differences between codon pairs.
"""
struct CodonTables
    aa::Vector{Char}
    ssites::Vector{Float64}
    syn::Matrix{Float64}        # NaN when no valid pathway (or a stop codon involved)
    nonsyn::Matrix{Float64}
end

_codon_bases(c::Int) = ((c - 1) ÷ 16, ((c - 1) ÷ 4) % 4, (c - 1) % 4)
_codon_from(b) = 16 * b[1] + 4 * b[2] + b[3] + 1

codon_tables(code::GeneticCode) = CODON_TABLES[code.id]

function _build_codon_tables(code::GeneticCode)
    aa = collect(code.table)
    # Synonymous sites: for each position, the fraction of the 3 possible
    # single-nucleotide changes that are synonymous. Changes to stop codons
    # are counted as nonsynonymous sites (Nei & Gojobori 1986).
    ss = zeros(64)
    for c in 1:64
        aa[c] == '*' && continue
        b = collect(_codon_bases(c))
        for p in 1:3, alt in 0:3
            alt == b[p] && continue
            bb = copy(b); bb[p] = alt
            aa[_codon_from(bb)] == aa[c] && (ss[c] += 1 / 3)
        end
    end
    syn = fill(NaN, 64, 64); non = fill(NaN, 64, 64)
    for c1 in 1:64, c2 in 1:64
        (aa[c1] == '*' || aa[c2] == '*') && continue
        b1 = collect(_codon_bases(c1)); b2 = collect(_codon_bases(c2))
        diffpos = [p for p in 1:3 if b1[p] != b2[p]]
        if isempty(diffpos)
            syn[c1, c2] = 0.0; non[c1, c2] = 0.0
            continue
        end
        s_tot = 0.0; n_tot = 0.0; npath = 0
        for order in _permutations(diffpos)
            cur = copy(b1); s = 0; n = 0; valid = true
            for p in order
                nxt = copy(cur); nxt[p] = b2[p]
                a_cur = aa[_codon_from(cur)]; a_nxt = aa[_codon_from(nxt)]
                if a_nxt == '*'
                    valid = false; break
                end
                a_cur == a_nxt ? (s += 1) : (n += 1)
                cur = nxt
            end
            valid || continue
            s_tot += s; n_tot += n; npath += 1
        end
        if npath > 0
            syn[c1, c2] = s_tot / npath; non[c1, c2] = n_tot / npath
        end
    end
    CodonTables(aa, ss, syn, non)
end

function _permutations(v::Vector{Int})
    length(v) <= 1 && return [copy(v)]
    out = Vector{Vector{Int}}()
    for i in eachindex(v)
        rest = [v[j] for j in eachindex(v) if j != i]
        for p in _permutations(rest)
            push!(out, vcat(v[i], p))
        end
    end
    out
end

const CODON_TABLES = Dict{Int,CodonTables}(id => _build_codon_tables(c) for (id, c) in GENETIC_CODES)

"""
Accumulator for codon comparisons of one sequence against another (or a
consensus).
"""
mutable struct CodonCompare
    n_codons::Int
    diff_pos::Vector{Int}       # nucleotide differences at codon positions 1..3
    n_diff_codons::Int
    n_multi::Int
    syn::Float64
    nonsyn::Float64
    n_unclassified::Int         # stop codon involved, or no valid pathway
    s_sites::Float64
    n_sites::Float64
end
CodonCompare() = CodonCompare(0, zeros(Int, 3), 0, 0, 0.0, 0.0, 0, 0.0, 0.0)

@inline function add_codon_pair!(cc::CodonCompare, t::CodonTables, c1::Int, c2::Int)
    cc.n_codons += 1
    b1 = _codon_bases(c1); b2 = _codon_bases(c2)
    nd = 0
    for p in 1:3
        if b1[p] != b2[p]
            cc.diff_pos[p] += 1; nd += 1
        end
    end
    if t.aa[c1] != '*' && t.aa[c2] != '*'
        cc.s_sites += (t.ssites[c1] + t.ssites[c2]) / 2
        cc.n_sites += 3 - (t.ssites[c1] + t.ssites[c2]) / 2
    end
    nd == 0 && return cc
    cc.n_diff_codons += 1
    nd > 1 && (cc.n_multi += 1)
    s = t.syn[c1, c2]
    if isnan(s)
        cc.n_unclassified += 1
    else
        cc.syn += s; cc.nonsyn += t.nonsyn[c1, c2]
    end
    cc
end

_jc(p) = p < 0.75 ? -0.75 * log(1 - 4p / 3) : NaN

"""
Summarise a `CodonCompare`. dN/dS (Nei-Gojobori with Jukes-Cantor
correction) is estimated only when `dnds = true` and its assumptions hold;
otherwise the values are `missing` and `dnds_reason` says why.
"""
function codon_summary(cc::CodonCompare; dnds::Bool = false, min_codons::Int = 30)
    nsn = cc.syn + cc.nonsyn
    base = (n_codons_compared = cc.n_codons, diff_pos1 = cc.diff_pos[1],
            diff_pos2 = cc.diff_pos[2], diff_pos3 = cc.diff_pos[3],
            diff_pos12 = cc.diff_pos[1] + cc.diff_pos[2],
            n_diff_codons = cc.n_diff_codons, n_multi_diff_codons = cc.n_multi,
            syn_diffs = cc.syn, nonsyn_diffs = cc.nonsyn,
            prop_nonsyn = nsn > 0 ? cc.nonsyn / nsn : missing,
            n_unclassified_codons = cc.n_unclassified,
            syn_sites = cc.s_sites, nonsyn_sites = cc.n_sites)
    dnds || return merge(base, (pN = missing, pS = missing, dN = missing, dS = missing,
                                dN_dS = missing, dnds_reason = "not requested"))
    pS = cc.s_sites > 0 ? cc.syn / cc.s_sites : NaN
    pN = cc.n_sites > 0 ? cc.nonsyn / cc.n_sites : NaN
    reason = if cc.n_codons < min_codons
        "fewer than $min_codons comparable codons"
    elseif isnan(pS) || isnan(pN)
        "no synonymous or nonsynonymous sites"
    elseif pS >= 0.75 || pN >= 0.75
        "p-distance ≥ 0.75: Jukes-Cantor correction undefined (saturation)"
    elseif cc.syn == 0
        "no synonymous differences: dN/dS not estimable"
    else
        ""
    end
    if reason != ""
        return merge(base, (pN = isnan(pN) ? missing : pN, pS = isnan(pS) ? missing : pS,
                            dN = missing, dS = missing, dN_dS = missing, dnds_reason = reason))
    end
    dN = _jc(pN); dS = _jc(pS)
    merge(base, (pN = pN, pS = pS, dN = dN, dS = dS, dN_dS = dN / dS,
                 dnds_reason = "estimated; between closely related sequences dN/dS " *
                               "reflects polymorphism, not selection (Kryazhimskiy & Plotkin 2008)"))
end

"""
    compare_codons(seq1, seq2; genetic_code = 5, frame = 1, dnds = false) -> NamedTuple

Codon-aware comparison of two **aligned** coding sequences (same length,
same frame). Codons with ambiguous bases or gaps are skipped
(`n_codons_skipped`).

Differences are counted at each codon position (`diff_pos1..3`).
Synonymous and nonsynonymous differences follow Nei & Gojobori (1986):
a codon differing at one position is classified directly; codons differing
at two or three positions are averaged over all mutational pathways that
avoid stop codons (fractional counts). Codon pairs involving a stop codon,
or with no valid pathway, are counted in `n_unclassified_codons` and not
classified. Differences at third positions are therefore **not** assumed
synonymous.

With `dnds = true`, pN/pS, Jukes-Cantor dN/dS and their ratio are
returned when estimable (≥ 30 codons, at least one synonymous difference,
p < 0.75); otherwise `missing` with `dnds_reason`. dN/dS is never the
raw ratio of difference counts.
"""
function compare_codons(seq1, seq2; genetic_code = 5, frame::Integer = 1, dnds::Bool = false)
    p1 = parse_sequence(seq1); p2 = parse_sequence(seq2)
    (p1 === missing || p2 === missing) && throw(ArgumentError("sequence is missing"))
    length(p1.mask) == length(p2.mask) || throw(ArgumentError(
        "sequences differ in length; codons can only be compared between aligned sequences"))
    t = codon_tables(MitoASV.genetic_code(genetic_code))
    cc = CodonCompare(); skipped = 0
    lead = frame - 1
    for k in 1:((length(p1.mask) - lead) ÷ 3)
        o = lead + 3(k - 1)
        c1 = codon_index(p1.mask[o + 1], p1.mask[o + 2], p1.mask[o + 3])
        c2 = codon_index(p2.mask[o + 1], p2.mask[o + 2], p2.mask[o + 3])
        if c1 == 0 || c2 == 0
            skipped += 1
            continue
        end
        add_codon_pair!(cc, t, c1, c2)
    end
    merge(codon_summary(cc; dnds), (n_codons_skipped = skipped,))
end

"""
    codon_usage(sequences; genetic_code = 5, frame = 1) -> DataFrame

Codon counts over all sequences (each sequence counted once; ambiguous and
gap codons skipped) and relative synonymous codon usage (RSCU, Sharp & Li
1987): the observed count divided by the mean count of the synonymous
codons. RSCU is `missing` for stop codons, single-codon amino acids and
amino acids that were never observed.
"""
function codon_usage(sequences; genetic_code = 5, frame::Integer = 1)
    code = MitoASV.genetic_code(genetic_code)
    counts = zeros(Int, 64)
    for s in sequences
        p = parse_sequence(s)
        p === missing && continue
        _count_codons!(counts, p.mask, frame)
    end
    _usage_table(counts, code)
end

function _count_codons!(counts, mask, frame)
    lead = frame - 1
    for k in 1:((length(mask) - lead) ÷ 3)
        o = lead + 3(k - 1)
        c = codon_index(mask[o + 1], mask[o + 2], mask[o + 3])
        c == 0 || (counts[c] += 1)
    end
    counts
end

function _usage_table(counts, code::GeneticCode)
    aa = collect(code.table)
    rscu = Vector{Union{Missing,Float64}}(missing, 64)
    for c in 1:64
        aa[c] == '*' && continue
        syn = findall(==(aa[c]), aa)
        tot = sum(counts[syn])
        (length(syn) == 1 || tot == 0) && continue
        rscu[c] = counts[c] / (tot / length(syn))
    end
    DataFrame(codon = [codon_string(c) for c in 1:64], amino_acid = aa, count = counts,
              n_synonymous_codons = [count(==(aa[c]), aa) for c in 1:64], rscu = rscu)
end
