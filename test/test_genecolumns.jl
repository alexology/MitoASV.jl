@testset "One column per gene and length check (SIMULATED data)" begin
    @test modal_length([300, 300, 303]) == 300
    @test modal_length([300, 303]) === nothing          # tie: no single expected length
    @test modal_length(Int[]) === nothing

    sim = simulate_asv_dataset(seed = 3)
    df = select(sim.data, Not(:sequence))                 # the usual input: gene columns only
    allowmissing!(df, [:ND2, :COI])
    i = findfirst(==("SIM_A01"), df.ASV)
    ref = df.COI[i]
    # in-frame internal deletion (3 nt), terminal truncation (3 nt), missing ND2
    push!(df, merge(NamedTuple(df[i, :]), (ASV = "SIM_DEL3", COI = ref[1:600] * ref[604:end])))
    push!(df, merge(NamedTuple(df[i, :]), (ASV = "SIM_TRUNC3", COI = ref[1:end - 3])))
    push!(df, merge(NamedTuple(df[i, :]), (ASV = "SIM_NOND2", ND2 = missing,
                                           COI = synonymous_change(ref, syn_site(ref, 40)))))
    # species beta has its own, shorter COI: its expected length must come from beta only
    beta = df.species .== "Simulatus beta"
    df.COI[beta] = [s[1:end - 6] for s in df.COI[beta]]

    res = assess_asvs(df; sequence_col = nothing, id_col = :ASV, genes = Dict(:ND2 => :ND2, :COI => :COI),
                      frame = 1, species_col = :species, abundance_col = :abundance)
    fl = res.flags
    has(id, f) = any((fl.id .== id) .& (fl.flag .== f))
    gfl(id, f, gene) = any((fl.id .== id) .& (fl.flag .== f) .& isequal.(fl.gene, gene))

    # groups are built separately for each gene
    @test all(occursin.("grouped separately for each gene", res.clusters.cluster_basis[.!ismissing.(res.clusters.cluster)]))
    @test Set(res.clusters.gene) == Set([:ND2, :COI])

    # length test
    ln = res.lengths
    row(id, gene) = findfirst((ln.id .== id) .& (ln.gene .== gene))
    @test ln.modal_length_nt[row("SIM_A02", :COI)] == length(ref)
    @test ln.modal_length_nt[row("SIM_C01", :COI)] == length(ref) - 6     # within species beta
    @test !any(has(id, :length_outlier) for id in df.ASV[beta])
    @test ln.length_difference[row("SIM_DEL3", :COI)] == -3
    @test ln.indel_location[row("SIM_DEL3", :COI)] == :internal
    @test gfl("SIM_DEL3", :length_outlier, :COI) && !gfl("SIM_DEL3", :invalid_frame, :COI)
    @test ln.indel_location[row("SIM_TRUNC3", :COI)] == :terminal
    @test gfl("SIM_TRUNC3", :length_outlier, :COI)
    @test gfl("SIM_FS1", :length_outlier, :COI) && gfl("SIM_FS1", :invalid_frame, :COI)
    @test ln.indel_location[row("SIM_FS1", :COI)] == :internal
    @test !has("SIM_A03", :length_outlier)
    @test ln.status[row("SIM_NOND2", :ND2)] == :not_evaluated
    @test any((res.tests.test .== :length) .& (res.tests.status .== :evaluated))

    # the other artefacts are still found with gene columns
    @test gfl("SIM_STOP1", :internal_stop, :COI)
    @test gfl("SIM_RAD1", :grantham_outlier, :COI)
    @test has("SIM_CHIM1", :gene_discordance)
    # a row without ND2 is still compared for COI
    aa = res.amino_acid_changes
    k = findfirst((aa.id .== "SIM_NOND2") .& (aa.gene .== :COI))
    @test aa.n_compared[k] !== missing && aa.n_compared[k] > 0
    @test res.gene_concordance.reason[findfirst(==("SIM_NOND2"), res.gene_concordance.id)] == "both genes are needed"

    # too few comparators in a species: length not evaluated, never passed
    small = assess_asvs(df[df.species .== "Simulatus beta", :][1:3, :]; sequence_col = nothing,
                        genes = [:COI => :COI], frame = 1, species_col = :species)
    @test all(small.lengths.status .== :not_evaluated)
    @test all(small.summary.classification .== "Insufficient evidence")

    # wide output carries per-gene length and cluster columns
    w = assess_asvs_df(df; sequence_col = nothing, id_col = :ASV, genes = Dict(:ND2 => :ND2, :COI => :COI),
                       frame = 1, species_col = :species)
    @test hasproperty(w, :COI_modal_length_nt) && hasproperty(w, :ND2_cluster)
    @test w.flag_length_outlier[findfirst(==("SIM_DEL3"), w.ASV)] === true
    @test_throws ArgumentError assess_asvs(df; sequence_col = nothing, genes = Dict(:COI => :nope), frame = 1)
end
