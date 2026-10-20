# Single-gene (COI) assessment with an unknown reading frame, on SIMULATED data.
# Run from the package folder:  julia --project=. examples/single_gene.jl
using MitoASV, DataFrames

sim = simulate_asv_dataset(seed = 2)
coi = sim.genes[:COI]
df = DataFrame(ASV = sim.data.ASV, species = sim.data.species,
               COI = [s[coi.start:min(coi.stop, end)] for s in sim.data.sequence])

# frame = nothing: the three frames are evaluated; a frame is used only if the
# dataset agrees on it (see results.parameters, row "frame_notes").
res = assess_asvs(df; sequence_col = :COI, id_col = :ASV, species_col = :species,
                  gene = :COI, frame = nothing)
println(res.parameters[res.parameters.parameter .== "frames_used", :value][1])
show(select(res.validation, :id, :frame, :frame_source, :n_internal_stops, :trailing_bases,
            :open_frames); allrows = true); println("\n")

# Stand-alone functions
p = translate_sequence(df.COI[1]; frame = 1)
println("mean KD of the first ASV: ", hydrophobicity_metrics(p).mean_kd)
rad = df.COI[findfirst(==("SIM_RAD1"), df.ASV)]          # same species as SIM_A01
println(compare_proteins(p, translate_sequence(rad; frame = 1)).n_substitutions,
        " substitution(s) between the first ASV and the radical-substitution ASV")

# One wide row per ASV
wide = assess_asvs_df(df; sequence_col = :COI, id_col = :ASV, species_col = :species,
                      gene = :COI, frame = 1)
show(select(wide, :ASV, :classification, :flag_internal_stop, :flag_grantham_outlier,
            :COI_mean_kd); allrows = true); println()
