@testset "DataFrame interface (single gene)" begin
    orf = simulated_orf(80; seed = 11)
    seqs = Any[orf, setbase(orf, 30, 'C'), setcodon(orf, 10, "TAA"), missing, "",
               setbase(orf, 50, '?'), orf * "A"]
    df = DataFrame(ASV = ["a", "b", "c", "d", "e", "f", "g"], seq = seqs, extra = 1:7)
    before = deepcopy(df)
    res = assess_asvs(df; sequence_col = :seq, id_col = :ASV, gene = :COI, frame = 1)
    @test isequal(df, before)                        # input untouched
    @test nrow(res.summary) == 7
    @test res.summary.id == df.ASV
    v = res.validation
    @test v.gene == fill(:COI, 7)
    @test v.n_internal_stops[3] == 1
    @test v.status[4] == :missing_sequence
    @test v.status[5] == :empty_sequence
    @test v.n_invalid_symbols[6] == 1
    @test v.trailing_bases[7] == 1
    @test v.n_internal_stops[7] == 0
    fl = res.flags
    @test any((fl.id .== "c") .& (fl.flag .== :internal_stop))
    @test any((fl.id .== "d") .& (fl.flag .== :missing_sequence))
    @test any((fl.id .== "f") .& (fl.flag .== :invalid_symbols))
    @test !any((fl.id .== "g") .& (fl.flag .== :invalid_frame))   # terminal extra base ≠ frameshift
    @test res.lengths.indel_location[7] == :terminal
    @test res.summary.classification[4] == "Insufficient evidence"
    @test res.summary.classification[3] != "No detected anomaly"
    @test nrow(res.hydrophobicity) == 7
    @test res.hydrophobicity.mean_kd[4] === missing
    @test "genetic_code" in res.parameters.parameter
end

@testset "Frame consensus and unresolved frames" begin
    orf = simulated_orf(120; seed = 5)
    df = DataFrame(sequence = [orf, setbase(orf, 9, 'A'), setbase(orf, 99, 'G'),
                               setbase(orf, 150, 'T'), setcodon(orf, 40, "TAG")])
    res = assess_asvs(df; gene = :COI)
    # long ORFs: frame 1 is the only open frame of the intact sequences
    @test all(res.validation.frame_source .== :consensus)
    @test all(res.validation.frame .== 1)
    @test res.validation.n_internal_stops[5] == 1    # not rescued by another frame
    # short sequences with several open frames: no consensus, frame unresolved
    df2 = DataFrame(sequence = ["ATGATGATGATG", "ATGATGATGATA", "ATGATAATGATG"])
    r2 = assess_asvs(df2)
    @test all(r2.validation.frame_source .== :unresolved)
    @test all(ismissing, r2.validation.n_internal_stops)
    @test all(r2.summary.classification .== "Insufficient evidence")
    @test any(r2.flags.flag .== :frame_unresolved)
end

@testset "BioSequences input" begin
    using BioSequences
    orf = simulated_orf(30; seed = 2)
    df = DataFrame(s = [LongDNA{4}(orf), LongRNA{4}(replace(orf, 'T' => 'U'))])
    r = assess_asvs(df; sequence_col = :s, frame = 1)
    @test r.validation.protein[1] == r.validation.protein[2] == translate_sequence(orf)
end
