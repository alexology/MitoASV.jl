# Kyte-Doolittle hydrophobicity metrics.
#
# Hydrophobicity is a weak diagnostic feature: mitochondrial membrane
# proteins (ND2, COI) are hydrophobic, but a plausible profile does not prove
# that a protein is functional, and no universal threshold separates genuine
# from artifactual sequences.

"""
    hydrophobicity_profile(protein; window = 19) -> Vector{Union{Missing,Float64}}

Sliding-window mean Kyte-Doolittle hydropathy. Element `k` is the mean of
residues `k:(k + window - 1)`; windows that contain a residue without a
hydropathy value (`X`, `*`, `-`, `B`, `Z`, ...) are `missing`, so stops and
gaps are never averaged over. The result is empty when the protein is
shorter than `window`.
"""
function hydrophobicity_profile(protein::AbstractString; window::Integer = 19)
    window >= 1 || throw(ArgumentError("window must be ≥ 1"))
    v = [kd_value(c) for c in protein]
    n = length(v) - window + 1
    out = Vector{Union{Missing,Float64}}(missing, max(n, 0))
    n <= 0 && return out
    s = 0.0; bad = 0
    for i in 1:window
        isnan(v[i]) ? (bad += 1) : (s += v[i])
    end
    for k in 1:n
        if k > 1
            o = v[k - 1]; isnan(o) ? (bad -= 1) : (s -= o)
            a = v[k + window - 1]; isnan(a) ? (bad += 1) : (s += a)
        end
        bad == 0 && (out[k] = s / window)
    end
    out
end

"""
    hydrophobicity_metrics(protein; window = 19, window_threshold = 1.6,
                           hydrophobic_cutoff = 0.0) -> NamedTuple

Kyte-Doolittle summary of one translated sequence. Residue statistics use
the 20 standard residues only; `n_residues_scored` gives their number.

- `mean_kd`, `median_kd`, `sd_kd`: residue hydropathy statistics.
- `frac_hydrophobic`: fraction of residues with hydropathy > `hydrophobic_cutoff`.
- `n_windows`, `frac_windows_above`: number of complete windows and the
  fraction whose mean exceeds `window_threshold`. The default 1.6 with a
  19-residue window is the transmembrane-segment heuristic of Kyte &
  Doolittle (1982); it is descriptive, not a validity threshold.
- `max_window_kd`, `min_window_kd`.
- `reason`: why some values are `missing` (`""` if none are).
"""
function hydrophobicity_metrics(protein::AbstractString; window::Integer = 19,
                                window_threshold::Real = 1.6, hydrophobic_cutoff::Real = 0.0)
    vals = Float64[]
    for c in protein
        x = kd_value(c)
        isnan(x) || push!(vals, x)
    end
    n = length(vals)
    reasons = String[]
    n == 0 && push!(reasons, "no scorable residues")
    n == 1 && push!(reasons, "sd needs ≥ 2 residues")
    prof = collect(skipmissing(hydrophobicity_profile(protein; window)))
    nw = length(prof)
    nw == 0 && push!(reasons, "no complete $(window)-residue window without stops, gaps or X")
    (n_residues_scored = n,
     mean_kd = n > 0 ? mean(vals) : missing,
     median_kd = n > 0 ? median(vals) : missing,
     sd_kd = n > 1 ? std(vals) : missing,
     frac_hydrophobic = n > 0 ? count(>(hydrophobic_cutoff), vals) / n : missing,
     n_windows = nw,
     frac_windows_above = nw > 0 ? count(>(window_threshold), prof) / nw : missing,
     max_window_kd = nw > 0 ? maximum(prof) : missing,
     min_window_kd = nw > 0 ? minimum(prof) : missing,
     reason = join(reasons, "; "))
end

"""
    hydrophobicity_table(proteins; ids = 1:length(proteins), kwargs...) -> DataFrame

[`hydrophobicity_metrics`](@ref) for several proteins, one row each, for
comparisons among ASVs. `missing` proteins give a row of `missing` values.
"""
function hydrophobicity_table(proteins::AbstractVector; ids = 1:length(proteins), kwargs...)
    rows = map(proteins) do p
        p === missing ? _missing_hydro("missing protein") : hydrophobicity_metrics(p; kwargs...)
    end
    df = DataFrame(rows)
    insertcols!(df, 1, :id => collect(ids))
    df
end

# Metrics of one unique sequence as stored by the pipeline (`GeneData.hydro`):
# a concrete type for both the computed and the missing case.
const HydroMetrics = @NamedTuple{n_residues_scored::Union{Missing,Int}, mean_kd::Union{Missing,Float64},
                                 median_kd::Union{Missing,Float64}, sd_kd::Union{Missing,Float64},
                                 frac_hydrophobic::Union{Missing,Float64}, n_windows::Union{Missing,Int},
                                 frac_windows_above::Union{Missing,Float64}, max_window_kd::Union{Missing,Float64},
                                 min_window_kd::Union{Missing,Float64}, reason::String}

_missing_hydro(reason) = (n_residues_scored = missing, mean_kd = missing, median_kd = missing,
                          sd_kd = missing, frac_hydrophobic = missing, n_windows = missing,
                          frac_windows_above = missing, max_window_kd = missing,
                          min_window_kd = missing, reason = reason)
