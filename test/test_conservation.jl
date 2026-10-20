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
end
