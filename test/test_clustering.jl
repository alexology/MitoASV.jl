@testset "Sequence clustering" begin
    base = simulated_orf(100; seed = 31)          # 300 nt
    mut(s, ps) = foldl((x, p) -> setbase(x, p, x[p] == 'A' ? 'C' : 'A'), ps; init = s)
    a = base
    b = mut(a, [10, 50, 90, 130, 170, 210, 250, 290])           # 8/300 = 2.7% from a
    c = mut(b, [20, 60, 100, 140, 180, 220, 260, 280])          # 2.7% from b, 5.3% from a
    d = mut(a, collect(5:5:150))                                 # 10% from a
    cl = cluster_asvs([a, b, c, d]; threshold = 0.03)
    # connected components chain a-b-c
    @test cl.cluster[1] == cl.cluster[2] == cl.cluster[3]
    @test cl.cluster[4] != cl.cluster[1]
    @test cl.cluster_size == [3, 3, 3, 1]
    @test cl.nn_distance[1] ≈ 8 / 300
    @test cl.nn_distance[4] === missing
    # complete and average linkage do not chain
    for lk in (:complete, :average)
        h = cluster_asvs([a, b, c, d]; method = :hierarchical, linkage = lk, threshold = 0.03)
        @test length(unique(h.cluster)) == 3
        @test h.cluster[1] == h.cluster[2] || h.cluster[2] == h.cluster[3]
    end
    hs = cluster_asvs([a, b, c, d]; method = :hierarchical, linkage = :single, threshold = 0.03)
    @test hs.cluster == cl.cluster
    # unequal length (1-nt deletion) still links; disabled alignment does not
    del = a[1:149] * a[151:end]
    u = cluster_asvs([a, del]; threshold = 0.03)
    @test u.cluster[1] == u.cluster[2]
    u2 = cluster_asvs([a, del]; threshold = 0.03, align_unequal = false)
    @test u2.cluster[1] != u2.cluster[2]
    # identical sequences and missing values
    m = cluster_asvs([a, missing, a, ""])
    @test m.cluster[1] == m.cluster[3]
    @test m.cluster[2] === missing && m.cluster[4] === missing
    # protein level: synonymous variants are identical proteins
    s1 = synonymous_change(a, syn_site(a, 10))
    p = cluster_asvs([a, s1, d]; level = :protein, threshold = 0.0)
    @test p.cluster[1] == p.cluster[2]
    # threaded and serial runs agree
    seqs = [mut(a, collect(k:7:(k + 70))) for k in 1:40]
    @test cluster_asvs(seqs; threshold = 0.04).cluster ==
          cluster_asvs(seqs; threshold = 0.04, threaded = false).cluster
    # banded edit distance
    ed(x, y, k) = MitoASV.banded_edit_distance(parse_sequence(x).mask, parse_sequence(y).mask, k)
    @test ed("ACGTACGT", "ACGACGT", 3) == 1
    @test ed("AAAA", "TTTTTT", 2) == 3
    @test ed("ACNT", "ACGT", 1) == 0
    @test ed("GATTACA", "GCATGCU", 7) == 4
end
