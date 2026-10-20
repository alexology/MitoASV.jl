@testset "Comparisons within species (SIMULATED data)" begin
    sim = simulate_asv_dataset(seed = 3)
    df = sim.data
    @test Set(df.species) == Set(["Simulatus alpha", "Simulatus beta"])
    res = assess_asvs(df; id_col = :ASV, genes = sim.genes, species_col = :species,
                      abundance_col = :abundance)
    cl = res.clusters
    # no comparison group contains more than one species
    for gdf in groupby(dropmissing(cl, :cluster), :cluster)
        @test length(unique(gdf.species)) == 1
    end
    @test hasproperty(res.summary, :species) && res.summary.species == df.species
    @test all(r -> r.species == "Simulatus alpha" || r.species == "Simulatus beta",
              eachrow(res.group_statistics))
    # the simulated artefacts are still flagged within species alpha
    fl = res.flags
    has(id, f) = any((fl.id .== id) .& (fl.flag .== f))
    @test has("SIM_STOP1", :internal_stop)
    @test has("SIM_FS1", :invalid_frame)
    @test has("SIM_RAD1", :grantham_outlier)
    @test has("SIM_CHIM1", :gene_discordance)      # intraspecific chimera (lineages A × B)
    @test all(res.summary.classification[df.category .== "other_species"] .== "No detected anomaly")
    # concordance neighbours are always of the same species
    gc = res.gene_concordance
    @test all(gc.scope .== :species)
    sp_of = Dict(zip(df.ASV, df.species))
    for r in eachrow(gc)
        r.nn_gene1_id === missing && continue
        @test sp_of[r.nn_gene1_id] == r.species && sp_of[r.nn_gene2_id] == r.species
    end
    @test all(r -> sp_of[r.id_i] == sp_of[r.id_j] == r.species, eachrow(res.discordant_pairs))

    # clustering disabled: each species is exactly one group
    one = assess_asvs(df; genes = sim.genes, species_col = :species, clustering = false)
    @test length(unique(one.clusters.cluster)) == 2
    @test all(one.clusters.cluster_source .== :single_group)

    # the same sequence assigned to two species is analysed separately in each
    d2 = vcat(df, DataFrame(df[1, :]))
    d2.species[end] = "Simulatus beta"
    d2.ASV[end] = "SIM_A01_beta"
    r2 = assess_asvs(d2; id_col = :ASV, genes = sim.genes, species_col = :species)
    c2 = r2.clusters
    @test c2.cluster[1] != c2.cluster[end]
    @test c2.unit[1] != c2.unit[end]

    # missing species: no comparison group, comparative tests not evaluated
    d3 = allowmissing(copy(df), :species)
    d3.species[2] = missing
    r3 = assess_asvs(d3; id_col = :ASV, genes = sim.genes, species_col = :species)
    @test r3.clusters.cluster[2] === missing
    @test any((r3.flags.row .== 2) .& (r3.flags.flag .== :missing_species))
    @test r3.summary.classification[2] == "Insufficient evidence"
    @test occursin("species missing", r3.amino_acid_changes.reason[2])
    @test r3.gene_concordance.status[2] == :not_evaluated

    # user groups are combined with species
    d4 = copy(df); d4.grp = fill("g1", nrow(d4))
    r4 = assess_asvs(d4; genes = sim.genes, species_col = :species, cluster_col = :grp)
    @test length(unique(r4.clusters.cluster)) == 2

    # interspecific chimera (ND2 of alpha + COI of beta, labelled alpha): only a
    # dataset-wide concordance can find its second parent
    cs = sim.genes[:COI].start
    alpha1 = df.sequence[findfirst(==("SIM_A01"), df.ASV)]
    beta1 = df.sequence[findfirst(==("SIM_C01"), df.ASV)]
    d5 = vcat(df, DataFrame(df[1, :]))
    d5.ASV[end] = "SIM_XCHIM"; d5.sequence[end] = alpha1[1:cs - 1] * beta1[cs:end]
    all5 = assess_asvs(d5; id_col = :ASV, genes = sim.genes, species_col = :species,
                       concordance_scope = :all)
    @test any((all5.flags.id .== "SIM_XCHIM") .& (all5.flags.flag .== :gene_discordance))
    @test all(all5.gene_concordance.scope .== :all)
    @test_throws ArgumentError assess_asvs(df; genes = sim.genes, concordance_scope = :genus)
    @test_throws ArgumentError assess_asvs(df; genes = sim.genes, species_col = :taxon)
end
