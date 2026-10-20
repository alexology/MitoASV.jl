@testset "Grantham, BLOSUM62 and radical substitutions" begin
    @test grantham_distance('S', 'R') == 110
    @test grantham_distance('L', 'I') == 5
    @test grantham_distance('C', 'W') == 215
    @test grantham_distance('A', 'A') == 0
    @test grantham_distance('A', 'X') === missing
    @test grantham_distance('A', '*') === missing
    @test grantham_distance("L", "I") == 5
    @test blosum62_score('W', 'W') == 11
    @test blosum62_score('A', 'R') == -1
    @test blosum62_score('A', '-') === missing
    @test is_radical('C', 'W') === true
    @test is_radical('L', 'I') === false
    @test is_radical('L', 'X') === missing
    rc = RadicalCriteria(grantham_min = nothing, abs_delta_kd_min = 5)
    @test is_radical('I', 'D', rc)
    @test !is_radical('I', 'L', rc)
    @test_throws ArgumentError RadicalCriteria(grantham_min = nothing)
    rcall = RadicalCriteria(grantham_min = 100, blosum62_max = -3, combine = :all)
    @test is_radical('C', 'W', rcall) == (215 >= 100 && blosum62_score('C', 'W') <= -3)
end

@testset "compare_proteins" begin
    c = compare_proteins("MLKW", "MIKC")
    @test c.n_compared == 4
    @test c.n_substitutions == 2
    @test c.prop_different ≈ 0.5
    @test c.grantham_sum == 5 + 215
    @test c.grantham_max == 215
    @test c.grantham_mean ≈ 110
    @test c.n_radical == 1
    @test c.delta_kd_sum ≈ (4.5 - 3.8) + (2.5 - -0.9)
    @test c.blosum62_min == blosum62_score('W', 'C')
    @test nrow(c.substitutions) == 2
    @test c.substitutions.position == [2, 4]
    # gaps, stops and ambiguous residues are counted, not scored
    g = compare_proteins("M-K*XL", "MAKLAL")
    @test g.n_compared == 3
    @test g.n_gap_positions == 1
    @test g.n_stop_differences == 1
    @test g.n_uncomparable == 1
    @test g.n_substitutions == 0
    @test g.grantham_sum == 0
    @test g.grantham_mean === missing
    # identical
    i = compare_proteins("MLKW", "MLKW")
    @test i.n_substitutions == 0 && i.identity == 1.0 && !i.homology_warning
    # unaligned input is refused
    @test_throws ArgumentError compare_proteins("MLK", "MLKW")
    # low identity: homology warning
    @test compare_proteins("AAAA", "WWWW").homology_warning
end
