@testset "Robust anomaly scores" begin
    v = [1.0, 1.1, 0.9, 1.0, 1.05, 0.95, 5.0]
    s = robust_anomaly_scores(v)
    @test s.score[7] > 10
    @test all(abs.(s.score[1:6]) .< 3.5)
    @test s.dispersion_method[1] == "mad"
    @test all(s.reliable)
    # zero MAD: falls back to mean absolute deviation, then the resolution floor
    z = robust_anomaly_scores([0.0, 0, 0, 0, 0, 0, 0, 3]; resolution = 1.0)
    @test z.score[8] ≈ 3.0                     # (3 - 0) / floor 1, not infinite
    @test z.dispersion_method[8] == "resolution_floor"
    @test z.score[1] ≈ (0 - 0) / 1.0
    mz = robust_anomaly_scores([0.0, 0, 0, 0, 0, 0, 2, 2])
    @test mz.dispersion_method[1] == "mean_absolute_deviation"
    @test isfinite(mz.score[1])
    # zero-variance group: identical values score 0; a deviant without any
    # dispersion gets no score instead of an extreme one
    zv = robust_anomaly_scores(fill(2.0, 6))
    @test all(zv.score .== 0)
    zd = robust_anomaly_scores([2.0, 2, 2, 2, 2, 7])
    @test zd.score[6] === missing
    @test occursin("zero dispersion", zd.reason[6])
    # missing values and too few comparators
    m = robust_anomaly_scores([1.0, missing, 1.2, 0.8])
    @test m.score[2] === missing
    @test all(ismissing, m.score)              # only 2 comparators each
    @test !any(m.reliable)
    # leave-one-out: the focal value does not shift its own median
    lo = robust_anomaly_scores([0.0, 0, 0, 0, 10]; resolution = 1, loo = true)
    @test lo.median[5] == 0
end

@testset "Spearman correlation" begin
    @test spearman([1, 2, 3, 4], [10, 20, 30, 40]) ≈ 1
    @test spearman([1, 2, 3, 4], [4, 3, 2, 1]) ≈ -1
    @test spearman([1, 1, 1], [1, 2, 3]) === missing
    @test spearman([1, 2], [1, 2]) === missing
    @test spearman([1.0, 2, NaN, 4, 5], [1.0, 2, 3, 4, 5]) ≈ 1
end
