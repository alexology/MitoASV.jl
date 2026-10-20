# Combined ND2–COI assessment of a SIMULATED dataset (not real sequences).
# The input has one column per gene (ND2, COI) and a species column.
# Run from the package folder:  julia -t auto --project=. examples/simulated_nd2_coi.jl
using MitoASV, DataFrames

sim = simulate_asv_dataset(seed = 1)
println(sim.description)
df = select(sim.data, :ASV, :species, :ND2, :COI, :abundance, :sample, :category)

# genes: gene => column. species_col: every comparison stays within a species
# (groups, consensus, expected gene length, anomaly scores, ND2–COI concordance).
res = assess_asvs(df; sequence_col = nothing, id_col = :ASV,
                  genes = [:ND2 => :ND2, :COI => :COI], frame = 1, genetic_code = 5,
                  species_col = :species, abundance_col = :abundance, sample_col = :sample)
show(stdout, MIME("text/plain"), res); println("\n")

# Descriptive classification next to the simulated category (known only because
# the data are simulated).
cmp = hcat(select(df, :ASV, :species, :category), select(res.summary, :classification, :anomaly_flags))
show(cmp; allrows = true, truncate = 80); println("\n")

# Gene lengths against the length expected from the other sequences of the species
show(filter(r -> r.gene == :COI, res.lengths)[:, [:id, :length_nt, :modal_length_nt,
     :length_difference, :indel_location, :status]]; allrows = true); println("\n")

# The evidence behind each anomaly flag
anom = filter(r -> r.flag_type == :anomaly, res.flags)
show(select(anom, :id, :gene, :flag, :detail); allrows = true, truncate = 120); println()
