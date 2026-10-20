@testset "Within-group comparisons in assess_asvs (SIMULATED data)" begin
    orf = simulated_orf(150; seed = 41)
    vars = [orf]
    for k in (10, 30, 50, 70, 90, 110, 130)
        push!(vars, synonymous_change(orf, syn_site(orf, k)))
    end
    k = 75
    target = translate_codon(orf[3k - 2:3k]) == 'D' ? 'W' : 'D'
    radical = aa_change(orf, k, target)
    shifted = orf[1:200] * orf[202:end]                     # 1-nt deletion
    df = DataFrame(id = ["v$i" for i in 1:length(vars)], seq = vars)
    push!(df, ("radical", radical)); push!(df, ("shift", shifted))
    push!(df, ("dup", orf))                                 # identical to v1
    res = assess_asvs(df; sequence_col = :seq, id_col = :id, gene = :COI, frame = 1)
    @test all(res.clusters.cluster .== 1)
    aa = res.amino_acid_changes
    r = findfirst(==("radical"), aa.id)
    @test aa.n_substitutions[r] == 1
    @test aa.grantham_sum[r] == grantham_distance(translate_codon(orf[3k - 2:3k]), target)
    @test all(aa.n_substitutions[1:length(vars)] .== 0)
    cs = res.codon_statistics
    @test cs.syn_diffs[2] ≈ 1.0 && cs.nonsyn_diffs[2] ≈ 0.0
    @test cs.nonsyn_diffs[r] >= 1
    # duplicates are analysed once and get identical results
    @test isequal(aa[aa.id .== "dup", :n_compared], aa[aa.id .== "v1", :n_compared])
    # frameshift-like deletion: invalid_frame (relative to group modal length)
    fl = res.flags
    @test any((fl.id .== "shift") .& (fl.flag .== :invalid_frame))
    @test any((fl.id .== "shift") .& (fl.flag .== :not_comparable))
    @test cs.reason[findfirst(==("shift"), cs.id)] != ""
    g = res.group_statistics
    @test nrow(g) == 1 && g.status[1] == :ok
    @test g.var_frac_pos3[1] > g.var_frac_pos1[1]
    @test nrow(res.codon_usage) == 64
    # small groups: insufficient comparators, never "no anomaly"
    small = assess_asvs(df[1:3, :]; sequence_col = :seq, frame = 1)
    @test all(small.amino_acid_changes.n_compared .=== missing)
    @test any(small.flags.flag .== :insufficient_comparators)
    @test all(small.summary.classification .== "Insufficient evidence")
    # identical sequences only: one unique sequence -> insufficient
    ident = assess_asvs(DataFrame(sequence = fill(orf, 6)); frame = 1)
    @test all(ident.flags.flag .== :insufficient_comparators)
    # user-provided groups
    df.grp = [i <= 4 ? "A" : "B" for i in 1:nrow(df)]
    ug = assess_asvs(df; sequence_col = :seq, cluster_col = :grp, frame = 1)
    @test ug.clusters.cluster_label == df.grp
    @test all(ug.clusters.cluster_source .== :user)
    # clustering disabled: one comparison group
    one = assess_asvs(df; sequence_col = :seq, frame = 1, clustering = false)
    @test all(one.clusters.cluster .== 1)
end
