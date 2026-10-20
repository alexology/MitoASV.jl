# Pairwise distances between nucleotide (mask) and protein sequences.
#
# Conventions: `NaN` = not comparable (no shared informative position, or
# unequal lengths without alignment); `Inf` = above the requested cutoff
# (the exact value was not computed).

"""
p-distance between two equal-length mask vectors with pairwise deletion:
only positions where both bases are unambiguous are compared. Stops early
(returning `Inf`) once the distance must exceed `cutoff`.
"""
function hamming_pdistance(a::AbstractVector{UInt8}, b::AbstractVector{UInt8}; cutoff::Float64 = Inf)
    n = length(a)
    maxdiff = isfinite(cutoff) ? floor(Int, cutoff * n + 1e-9) : typemax(Int)
    diff = 0; comp = 0
    # Branch-free blocks (vectorisable); the cutoff is checked between blocks.
    i = 1
    @inbounds while i <= n
        stop = min(i + 511, n)
        for k in i:stop
            x = a[k]; y = b[k]
            ok = (count_ones(x) == 1) & (x < 0x10) & (count_ones(y) == 1) & (y < 0x10)
            comp += ok
            diff += ok & (x != y)
        end
        diff > maxdiff && return Inf
        i = stop + 1
    end
    comp == 0 && return NaN
    d = diff / comp
    d > cutoff ? Inf : d
end

"""
Banded edit distance (Ukkonen 1985) between ungapped mask vectors: unit
costs, ambiguous bases match any overlapping base. Returns the distance if
it is ≤ `k`, otherwise `k + 1`.
"""
function banded_edit_distance(a::AbstractVector{UInt8}, b::AbstractVector{UInt8}, k::Int)
    la, lb = length(a), length(b)
    abs(la - lb) > k && return k + 1
    la == 0 && return lb
    lb == 0 && return la
    big = k + 1
    w = 2k + 1
    prev = fill(big, w + 2); cur = fill(big, w + 2)
    # band index: column j stored at j - i + k + 2 (offset so that 1 is a sentinel)
    for j in 0:min(k, lb)
        prev[j + k + 2] = j
    end
    for i in 1:la
        fill!(cur, big)
        jlo = max(0, i - k); jhi = min(lb, i + k)
        rowmin = big
        for j in jlo:jhi
            idx = j - i + k + 2
            v = big
            if j == 0
                v = i
            else
                sub = (a[i] & b[j]) != 0x00 ? 0 : 1
                v = min(v, prev[idx] + sub)                    # (i-1, j-1)
                idx > 1 && (v = min(v, cur[idx - 1] + 1))      # (i, j-1)
            end
            idx + 1 <= w + 2 && (v = min(v, prev[idx + 1] + 1)) # (i-1, j)
            v = min(v, big)
            cur[idx] = v
            rowmin = min(rowmin, v)
        end
        rowmin > k && return big
        prev, cur = cur, prev
    end
    min(prev[lb - la + k + 2], big)
end

"""
Edit distance if ≤ `kmax`, else `kmax + 1`. The band starts narrow and is
doubled only when needed, so similar sequences cost O(n · |Δlength|).
"""
function capped_edit_distance(a::AbstractVector{UInt8}, b::AbstractVector{UInt8}, kmax::Int)
    k = min(max(abs(length(a) - length(b)) + 8, 16), kmax)
    while true
        d = banded_edit_distance(a, b, k)
        (d <= k || k >= kmax) && return min(d, kmax + 1)
        k = min(2k, kmax)
    end
end

"""
    sequence_distance(a, b; cutoff = Inf, align_unequal = true) -> Float64

Nucleotide distance between two parsed sequences (mask vectors).
Equal lengths: p-distance with pairwise deletion of gaps and ambiguous
positions (positional homology is assumed; aligned input satisfies it).
Unequal lengths: edit distance of the ungapped sequences divided by the
longer length, computed only with `align_unequal = true` (otherwise `NaN`).
"""
function sequence_distance(a::AbstractVector{UInt8}, b::AbstractVector{UInt8};
                           cutoff::Float64 = Inf, align_unequal::Bool = true)
    length(a) == length(b) && return hamming_pdistance(a, b; cutoff)
    align_unequal || return NaN
    ua = ungapped(a); ub = ungapped(b)
    L = max(length(ua), length(ub))
    L == 0 && return NaN
    k = isfinite(cutoff) ? floor(Int, cutoff * L + 1e-9) : L
    ed = capped_edit_distance(ua, ub, k)
    ed > k ? Inf : ed / L
end

"""
p-distance between two equal-length proteins over positions where both
residues are standard amino acids. `NaN` if lengths differ or no position
is comparable.
"""
function protein_pdistance(a::AbstractString, b::AbstractString; cutoff::Float64 = Inf)
    ncodeunits(a) == ncodeunits(b) || return NaN
    ca = codeunits(a); cb = codeunits(b)
    diff = 0; comp = 0
    maxdiff = isfinite(cutoff) ? floor(Int, cutoff * length(ca) + 1e-9) : typemax(Int)
    @inbounds for i in eachindex(ca)
        x = ca[i]; y = cb[i]
        if aa_index(x) != 0 && aa_index(y) != 0
            comp += 1
            if x != y
                diff += 1
                diff > maxdiff && return Inf
            end
        end
    end
    comp == 0 && return NaN
    d = diff / comp
    d > cutoff ? Inf : d
end
