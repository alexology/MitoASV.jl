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

    # the global strand applies to genes without their own, as the frame does
    rc = DataFrame(COI = ["ATGAAACCC"], species = ["s"])
    kw = (sequence_col = nothing, frame = 1, species_col = :species)
    minus = assess_asvs(rc; genes = Dict(:COI => (column = :COI, strand = :-)), kw...).validation.protein
    @test minus == ["GFH"]
    @test assess_asvs(rc; genes = Dict(:COI => :COI), strand = :-, kw...).validation.protein == minus
    @test assess_asvs(rc; genes = Dict(:COI => (column = :COI,)), strand = :-, kw...).validation.protein == minus
    @test assess_asvs(rc; genes = Dict(:COI => (column = :COI, strand = :+)), strand = :-,
                      kw...).validation.protein == ["MKP"]
end

@testset "Protein clustering respects the selected gene" begin
    df = DataFrame(g1 = ["ATGATGATG", "ATGATGATG"],
                   g2 = ["GCTGCTGCT", "GGTGGTGGT"])
    kw = (sequence_col = nothing, genes = Dict(:G1 => :g1, :G2 => :g2),
          frame = 1, cluster_level = :protein, cluster_threshold = 0.0,
          concordance = false, threaded = false)
    selected = assess_asvs(df; kw..., cluster_on = :G1)
    @test all(==(1), selected.clusters.cluster)
    @test all(==("protein translation of gene G1"), selected.clusters.cluster_basis)
    other = assess_asvs(df; kw..., cluster_on = :G2)
    @test length(unique(other.clusters.cluster)) == 2
    @test_throws ArgumentError assess_asvs(df; kw..., cluster_on = :unknown)
    # A missing unselected gene must not prevent clustering on the chosen one.
    allowmissing!(df, :g2)
    df.g2[2] = missing
    @test all(==(1), assess_asvs(df; kw..., cluster_on = :G1).clusters.cluster)
end

@testset "Wide output preserves existing prefixed columns" begin
    df = DataFrame(sequence = ["ATGATGATG"], classification = ["original"],
                   mito_classification = ["provenance"],
                   mito_mito_classification = ["another source"])
    original = copy(df)
    result = assess_asvs_df(df; frame = 1, comparative = false,
                           concordance = false, threaded = false)
    @test df == original
    @test result[:, names(df)] == original
    @test result.mito_mito_mito_classification == ["Insufficient evidence"]
    @test length(unique(names(result))) == ncol(result)
end
