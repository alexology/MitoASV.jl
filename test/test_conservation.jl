@testset "Reference-free conservation" begin
    core = "MLKWVIFAGSTLLIVAGMW"
    edit(s, i, c) = String(setindex!(collect(s), c, i))
    others = [core, core, core, edit(core, 3, 'R'), edit(core, 8, 'V')]
    radical = edit(core, 6, 'D')          # I -> D at a site invariant in the others
    res = conservation_analysis(vcat(others, [radical]))
    @test nrow(res.positions) == length(core)
    @test res.positions.conservation[1] ≈ 1.0
    @test res.positions.consensus[6] == 'I'
    sc = res.sequences
    @test sc.n_substitutions[6] == 1
    @test sc.n_radical_invariant[6] == 1
    @test sc.n_unseen_residues[6] == 1
    @test sc.grantham_sum[6] == grantham_distance('I', 'D')
    @test sc.cw_grantham_sum[6] ≈ grantham_distance('I', 'D')     # conservation 1 without it
    # leave-one-out: residues shared with others are not unseen
    @test sc.n_unseen_residues[1] == 0
    # the conservative K->R variant differs from its leave-one-out consensus at one position
    @test sc.n_substitutions[4] == 1
    @test sc.n_radical_invariant[4] == 0
    # insufficient sequences
    small = conservation_analysis(others[1:3])
    @test isempty(small.positions)
    @test occursin("fewer than", small.reason)
    # missing entries are skipped, unaligned refused
    withmissing = conservation_analysis(vcat(others, [missing, radical]))
    @test withmissing.sequences.reason[6] == "missing sequence"
    @test withmissing.sequences.n_radical_invariant[7] == 1
    @test_throws ArgumentError conservation_analysis(["MLK", "MLKW", "MLK", "MLK", "MLK"])
    # Henikoff weights down-weight redundant sequences
    w = henikoff_weights([Vector{UInt8}(codeunits(s)) for s in vcat(others, [radical])])
    @test w[1] ≈ w[2] ≈ w[3]
    @test w[6] > w[1]
    @test sum(w) ≈ 6
    hw = conservation_analysis(vcat(others, [radical]); weights = :henikoff)
    @test hw.n_effective < 6
    uw = conservation_analysis(vcat(others, [radical]); weights = :unique)
    @test uw.weights[2] == 0.0
    @test_throws ArgumentError conservation_analysis(others; weights = [1, 2])

    # Rescaling every weight changes nothing (no false ties, no empty consensus)
    s = ["A", "A", "A", "A", "V"]
    ref = conservation_analysis(s; weights = ones(5)).sequences
    for scale in (1e-10, 1e-15, 1e15)
        r = conservation_analysis(s; weights = fill(scale, 5)).sequences
        @test isequal(r.n_compared, ref.n_compared)
        @test isequal(r.grantham_sum, ref.grantham_sum)
    end
    # `:unique` counts distinct sequences and scores identical rows alike
    r = conservation_analysis(s; weights = :unique, min_sequences = 5)
    @test r.reason == "fewer than 5 distinct sequences (2)"
    r = conservation_analysis(["AK", "AK", "AK", "VK", "VR", "AR"]; weights = :unique, min_sequences = 3)
    @test r.reason == ""
    @test r.positions.n_informative == [4, 4]
    @test r.weights == [1.0, 0.0, 0.0, 1.0, 1.0, 1.0]
    @test allequal(r.sequences.n_compared[1:3]) && allequal(r.sequences.grantham_sum[1:3])
    @test r.sequences.grantham_sum[1] == 90   # AK against the consensus VR of VK, VR, AR
end
