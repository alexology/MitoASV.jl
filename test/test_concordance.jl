@testset "Gene concordance (SIMULATED recombinant)" begin
    a1 = simulated_orf(150; seed = 51); a2 = simulated_orf(200; seed = 52)
    # two lineages differing at ~10% of codons (synonymous changes)
    div(s, step) = foldl((x, k) -> something(try_synonymous(x, k), x), 1:step:(length(s) ÷ 3); init = s)
    b1 = div(a1, 8); b2 = div(a2, 8)
    varsA = [(synonymous_change(a1, syn_site(a1, 10i)), synonymous_change(a2, syn_site(a2, 12i))) for i in 1:5]
    varsB = [(synonymous_change(b1, syn_site(b1, 11i)), synonymous_change(b2, syn_site(b2, 13i))) for i in 1:5]
    chim = (varsA[1][1], varsB[1][2])
    all1 = vcat(first.(varsA), first.(varsB), [chim[1]])
    all2 = vcat(last.(varsA), last.(varsB), [chim[2]])
    gc = gene_concordance(all1, all2; k = 3)
    per = gc.per_asv
    @test per.discordant[end]
    @test per.nn_gene1[end] in 1:5 && per.nn_gene2[end] in 6:10
    @test per.two_parent_improvement[end] > 0.01
    @test count(per.discordant[1:10]) == 0
    @test gc.global_spearman_nt > 0.3
    @test per.knn_jaccard[end] < per.knn_jaccard[1]
    @test nrow(gc.pairs) >= 1
    # the genes may evolve at different rates: scaling one gene's distances
    # (all ranks unchanged) does not change rank-based results
    @test gc.global_spearman_nt == gene_concordance(all1, all2; k = 3).global_spearman_nt
    tr = gene_concordance(all1, all2; compare_trees = true)
    @test tr.tree_spearman isa Float64
end

@testset "Combined ND2–COI analysis (SIMULATED dataset)" begin
    sim = MitoASV.simulate_asv_dataset(seed = 3)
    df = sim.data
    @test all(df.simulated)
    res = assess_asvs(df; sequence_col = :sequence, id_col = :ASV, genes = sim.genes,
                      abundance_col = :abundance, sample_col = :sample)
    sm = res.summary
    cls(id) = sm.classification[findfirst(==(id), sm.id)]
    fl = res.flags
    has(id, f) = any((fl.id .== id) .& (fl.flag .== f))
    @test has("SIM_STOP1", :internal_stop)
    @test has("SIM_FS1", :invalid_frame)
    @test has("SIM_RAD1", :grantham_outlier)
    @test has("SIM_CHIM1", :gene_discordance)
    @test !has("SIM_NS1", :grantham_outlier)
    for id in ("SIM_STOP1", "SIM_FS1", "SIM_RAD1", "SIM_CHIM1")
        @test cls(id) in ("One or more anomalies", "Multiple independent anomalies")
    end
    # plausible variants: no anomaly flags at all
    plaus = df.ASV[df.category .== "plausible_variant"]
    @test all(cls(id) == "No detected anomaly" for id in plaus)
    # the divergent lineage forms its own group and is not flagged
    @test all(!has(id, :hydrophobicity_outlier) && !has(id, :gene_discordance)
              for id in df.ASV[df.category .== "divergent_plausible"])
    cl = res.clusters[res.clusters.gene .== :COI, :]
    @test length(unique(skipmissing(cl.cluster[df.lineage .== "A"]))) >= 1
    @test cl.cluster[findfirst(==("SIM_B01"), df.ASV)] != cl.cluster[findfirst(==("SIM_A01"), df.ASV)]
    # gene-specific tables
    @test Set(res.validation.gene) == Set([:ND2, :COI])
    @test nrow(res.validation) == 2nrow(df)
    @test all(res.gene_concordance.status .== :evaluated)
    # anomalies record raw value, group median, dispersion, score and reliability
    an = res.anomalies
    @test all(c -> hasproperty(an, c), (:value, :group_median, :group_dispersion, :robust_score,
                                        :reliable, :threshold, :dispersion_method))
    @test all(isfinite, skipmissing(an.robust_score))
    # wide output: one row per input row, input columns preserved
    w = assess_asvs_df(df; sequence_col = :sequence, id_col = :ASV, genes = sim.genes)
    @test nrow(w) == nrow(df)
    @test w.ASV == df.ASV && w.sequence == df.sequence
    @test hasproperty(w, :COI_mean_kd) && hasproperty(w, :ND2_n_internal_stops)
    @test w.flag_internal_stop[findfirst(==("SIM_STOP1"), w.ASV)] === true
    @test w.flag_internal_stop[findfirst(==("SIM_A01"), w.ASV)] === false
    # gene-specific anomaly thresholds are recorded
    r2 = assess_asvs(df; sequence_col = :sequence, genes = sim.genes,
                     anomaly_threshold = Dict(:ND2 => 4.0, :COI => 3.0))
    @test occursin("ND2 => 4.0", r2.parameters.value[r2.parameters.parameter .== "anomaly_threshold"][1])
    @test all(r2.anomalies.threshold[(r2.anomalies.gene .== :ND2) .& (r2.anomalies.method .!= "rule")] .== 4.0)
    # pre-extracted gene columns
    df2 = DataFrame(nd2 = [s[sim.genes[:ND2].start:sim.genes[:ND2].stop] for s in df.sequence],
                    coi = [length(s) >= sim.genes[:COI].stop ? s[sim.genes[:COI].start:sim.genes[:COI].stop] :
                           s[sim.genes[:COI].start:end] for s in df.sequence])
    r3 = assess_asvs(df2; sequence_col = nothing,
                     genes = [:ND2 => (column = :nd2, frame = 1), :COI => (column = :coi, frame = 1)])
    @test r3.validation.n_internal_stops == res.validation.n_internal_stops
    # reverse strand gene
    rc = DataFrame(seq = [MitoASV.reverse_complement(s) for s in df2.coi])
    r4 = assess_asvs(rc; sequence_col = :seq, gene = :COI, frame = 1, strand = :-)
    @test r4.validation.protein == res.validation.protein[res.validation.gene .== :COI]
end

@testset "Determinism and disabled analyses" begin
    sim = MitoASV.simulate_asv_dataset(seed = 5, n_variants = 6)
    a = assess_asvs(sim.data; genes = sim.genes)
    b = assess_asvs(sim.data; genes = sim.genes, threaded = false)
    @test isequal(a.summary, b.summary)
    @test isequal(a.anomalies, b.anomalies)
    @test isequal(a.gene_concordance, b.gene_concordance)
    c = assess_asvs(sim.data; genes = sim.genes, comparative = false, concordance = false)
    @test nrow(c.anomalies) == 0
    @test all(c.summary.classification .!= "No detected anomaly")    # missing ≠ passed
    @test any(c.tests.reason .== "disabled (comparative = false)")
    @test all(c.gene_concordance.status .== :not_evaluated)
end
