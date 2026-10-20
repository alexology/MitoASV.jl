@testset "Pipeline edge cases (SIMULATED data)" begin
    orf = simulated_orf(150; seed = 61)
    syn = [orf; [synonymous_change(orf, syn_site(orf, k)) for k in (10, 30, 50, 70, 90, 110, 130)]]
    # zero-variance group: synonymous variants share one protein, so protein
    # metrics have zero spread; scores are 0, never extreme
    z = assess_asvs(DataFrame(sequence = syn); frame = 1)
    kd = z.anomalies[z.anomalies.metric .== :mean_kd, :]
    @test all(kd.robust_score .== 0)
    @test all(kd.dispersion_method .== "resolution_floor")
    @test all(z.summary.classification .== "No detected anomaly")
    # ambiguous bases inside a group member: codon skipped, still compared
    amb = setbase(syn[2], 61, 'N')
    a = assess_asvs(DataFrame(sequence = [syn; amb]); frame = 1)
    v = a.validation
    @test v.n_ambiguous_codons[end] == 1
    @test a.codon_statistics.n_codons_compared[end] == a.codon_statistics.n_codons_compared[1] - 1
    # hierarchical and protein-level grouping inside assess_asvs
    h = assess_asvs(DataFrame(sequence = syn); frame = 1, cluster_method = :hierarchical)
    @test all(h.clusters.cluster .== 1)
    @test all(h.clusters.cluster_source .== :hierarchical)
    p = assess_asvs(DataFrame(sequence = syn); frame = 1, cluster_level = :protein, cluster_threshold = 0.0)
    @test all(p.clusters.cluster .== 1)
    # multi-gene with a truncated sequence and a missing row
    sim = simulate_asv_dataset(seed = 7, n_variants = 6)
    df = copy(sim.data)
    short = df.sequence[1][1:(sim.genes[:COI].start - 10)]
    push!(df, (ASV = "SIM_SHORT", sequence = short, species = "Simulatus alpha", category = "truncated", lineage = "A",
               abundance = 3, sample = "S1", simulated = true, ND2 = short[1:sim.genes[:ND2].stop], COI = ""))
    allowmissing!(df, [:sequence, :ND2, :COI])
    push!(df, (ASV = "SIM_MISSING", sequence = missing, species = "Simulatus alpha", category = "missing", lineage = "A",
               abundance = 1, sample = "S2", simulated = true, ND2 = missing, COI = missing))
    r = assess_asvs(df; id_col = :ASV, genes = sim.genes)
    fl = r.flags
    @test any((fl.id .== "SIM_SHORT") .& (fl.flag .== :gene_absent) .& isequal.(fl.gene, :COI))
    @test any((fl.id .== "SIM_MISSING") .& (fl.flag .== :missing_sequence))
    sm = r.summary
    @test sm.classification[sm.id .== "SIM_SHORT"] == ["Insufficient evidence"]
    @test sm.classification[sm.id .== "SIM_MISSING"] == ["Insufficient evidence"]
    @test r.validation.status[(r.validation.id .== "SIM_SHORT") .& (r.validation.gene .== :ND2)] == [:ok]
    @test nrow(r.summary) == nrow(df)
    # long format: the same sequence in several samples/replicates is analysed once
    long = vcat(sim.data, sim.data)
    long.replicate = repeat([1, 2]; inner = nrow(sim.data))
    lr = assess_asvs(long; id_col = :ASV, genes = sim.genes, abundance_col = :abundance,
                     sample_col = :sample, replicate_col = :replicate)
    @test all(lr.summary.n_rows_same_sequence .>= 2)
    @test all(lr.summary.n_replicates_detected .>= 2)
    @test isequal(lr.summary.classification[1:nrow(sim.data)],
                  lr.summary.classification[(nrow(sim.data) + 1):end])
    # argument checks
    @test_throws ArgumentError assess_asvs(sim.data; sequence_col = :nope)
    @test_throws ArgumentError assess_asvs(sim.data; genes = sim.genes, cluster_method = :kmeans)
    @test_throws ArgumentError GeneSpec(:COI; start = 0)
    @test_throws ArgumentError GeneSpec(:COI; start = 10, stop = 5)
    @test_throws ArgumentError GeneSpec(:COI; frame = 4)
end
