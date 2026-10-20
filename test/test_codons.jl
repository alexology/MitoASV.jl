@testset "Synonymous and nonsynonymous substitutions (NG86, table 5)" begin
    syn = compare_codons("TTT", "TTC")
    @test syn.syn_diffs == 1 && syn.nonsyn_diffs == 0
    non = compare_codons("TTT", "TTA")
    @test non.syn_diffs == 0 && non.nonsyn_diffs == 1
    # a third-position difference that is NOT synonymous (Ile -> Met in table 5)
    third = compare_codons("ATT", "ATG")
    @test third.diff_pos3 == 1
    @test third.nonsyn_diffs == 1 && third.syn_diffs == 0
    # synonymous in table 5 (ATA = Met) but not in table 1 (ATA = Ile)
    @test compare_codons("ATG", "ATA").syn_diffs == 1
    @test compare_codons("ATG", "ATA"; genetic_code = 1).nonsyn_diffs == 1
    # two differences: average over both pathways (CTT -> TTA, both Leu)
    two = compare_codons("CTT", "TTA")
    @test two.n_multi_diff_codons == 1
    @test two.syn_diffs ≈ 1.0 && two.nonsyn_diffs ≈ 1.0
    # pathway through a stop codon is excluded (TAT -> TGA via TAA)
    stp = compare_codons("TAT", "TGA")
    @test stp.syn_diffs ≈ 0 && stp.nonsyn_diffs ≈ 2
    # stop codon involved: unclassified, not counted
    u = compare_codons("TAA", "TAT")
    @test u.n_unclassified_codons == 1
    @test u.syn_diffs == 0 && u.nonsyn_diffs == 0
    @test u.prop_nonsyn === missing
    # ambiguous codons are skipped; unaligned refused
    @test compare_codons("TTNTTT", "TTTTTC").n_codons_skipped == 1
    @test_throws ArgumentError compare_codons("TTT", "TTTA")
    # synonymous sites of TTT in table 5: only TTC is synonymous -> 1/3
    t = MitoASV.codon_tables(genetic_code(5))
    @test t.ssites[MitoASV.codon_index(0x08, 0x08, 0x08)] ≈ 1 / 3
end

@testset "dN/dS only when estimable" begin
    orf = simulated_orf(100; seed = 21)
    r = compare_codons(orf, orf; dnds = true)
    @test r.dN_dS === missing
    @test occursin("no synonymous", r.dnds_reason)
    @test compare_codons(orf, orf).dnds_reason == "not requested"
    v = orf
    for k in (5, 20, 40, 60)
        v = synonymous_change(v, syn_site(v, k))
    end
    v = aa_change(v, 80, translate_codon(v[238:240]) == 'W' ? 'C' : 'W')
    r2 = compare_codons(orf, v; dnds = true)
    @test r2.syn_diffs >= 4
    @test r2.dN_dS isa Float64 && r2.dN_dS > 0
    @test r2.dS ≈ -0.75 * log(1 - 4 / 3 * r2.pS)
    short = compare_codons(orf[1:30], v[1:30]; dnds = true)
    @test short.dN_dS === missing
end

@testset "Codon usage and RSCU" begin
    cu = codon_usage(["TTTTTTTTC", "TGGATG"])
    @test cu.count[cu.codon .== "TTT"] == [2]
    @test cu.rscu[cu.codon .== "TTT"][1] ≈ 2 / (3 / 2)
    @test cu.rscu[cu.codon .== "TAA"][1] === missing        # stop
    @test cu.rscu[cu.codon .== "GCT"][1] === missing        # Ala never observed
    @test cu.n_synonymous_codons[cu.codon .== "TGG"] == [2] # TGA/TGG = W in table 5
end
