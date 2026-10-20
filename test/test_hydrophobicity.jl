@testset "Kyte-Doolittle hydrophobicity" begin
    @test KYTE_DOOLITTLE['I'] == 4.5
    @test KYTE_DOOLITTLE['R'] == -4.5
    prof = hydrophobicity_profile("IIIRRR"; window = 3)
    @test length(prof) == 4
    @test prof[1] ≈ 4.5
    @test prof[2] ≈ (4.5 + 4.5 - 4.5) / 3
    @test prof[4] ≈ -4.5
    # windows over unscorable residues are missing, never averaged across a stop
    p2 = hydrophobicity_profile("II*II"; window = 2)
    @test isequal(p2, [4.5, missing, missing, 4.5])
    @test isempty(hydrophobicity_profile("II"; window = 19))
    m = hydrophobicity_metrics("IIIRRR"; window = 3, window_threshold = 1.6)
    @test m.n_residues_scored == 6
    @test m.mean_kd ≈ 0.0
    @test m.median_kd ≈ 0.0
    @test m.sd_kd ≈ std([4.5, 4.5, 4.5, -4.5, -4.5, -4.5])
    @test m.frac_hydrophobic ≈ 0.5
    @test m.n_windows == 4
    @test m.frac_windows_above ≈ 0.25          # windows III (4.5), IIR (1.5), IRR, RRR
    @test m.max_window_kd ≈ 4.5
    @test m.min_window_kd ≈ -4.5
    short = hydrophobicity_metrics("IL"; window = 19)
    @test short.mean_kd ≈ (4.5 + 3.8) / 2
    @test short.frac_windows_above === missing
    @test occursin("window", short.reason)
    empty = hydrophobicity_metrics("X*"; window = 3)
    @test empty.mean_kd === missing && empty.sd_kd === missing
    t = hydrophobicity_table(["IIIRRR", missing]; window = 3)
    @test nrow(t) == 2
    @test t.mean_kd[2] === missing
end
