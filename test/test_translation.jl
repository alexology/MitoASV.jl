@testset "Genetic code and translation" begin
    # Table 5 (invertebrate mitochondrial): TGA = W, AGA/AGG = S, ATA = M.
    s = "ATGTGAAGAAGGATATAA"
    @test translate_sequence(s; genetic_code = 5) == "MWSSM*"
    @test translate_sequence(s; genetic_code = 1) == "M*RRI*"
    @test translate_sequence(s; genetic_code = 2) == "MW**M*"
    @test translate_codon("TGA") == 'W'
    @test translate_codon("AGA"; genetic_code = 1) == 'R'
    @test genetic_code(5).name == "Invertebrate Mitochondrial"
    @test_throws ArgumentError genetic_code(99)
    # all 64 codons of table 5
    tbl = genetic_code(5).table
    for (k, c) in enumerate(a * b * d for a in "TCAG" for b in "TCAG" for d in "TCAG")
        @test translate_codon(c) == tbl[k]
    end
    # lowercase, RNA and frames
    @test translate_sequence("augugaaga") == "MWS"
    @test translate_sequence("CATGTGA"; frame = 2) == "MW"
    @test translate_sequence("CCATGTGA"; frame = 3) == "MW"
    @test translate_sequence(MitoASV.reverse_complement("ATGTGAAGA"); strand = :-) == "MWS"
    @test_throws ArgumentError translate_sequence("ATG"; frame = 4)
    @test translate_sequence(missing) === missing
end

@testset "Ambiguous nucleotides" begin
    # CTN is always Leu; TGR = TGA/TGG = W in table 5 but * / W in table 1
    @test translate_sequence("CTN") == "L"
    @test translate_sequence("TGR") == "W"
    @test translate_sequence("TGR"; genetic_code = 1) == "X"
    @test translate_sequence("NNN") == "X"
    r = coding_validation("ATGCTNTGRTAR"; frame = 1)
    @test r.record.n_ambiguous_codons == 3
    @test r.record.n_unresolved_codons == 0     # TAR = TAA/TAG = stop in table 5
    @test r.record.terminal_stop
    r1 = coding_validation("ATGTRAGGG"; frame = 1)   # TRA = TAA/TGA: stop or W
    @test r1.protein == "MXG"
    @test r1.record.n_possible_stops == 1
    @test r1.record.n_unresolved_codons == 1
    @test r1.n_internal_stops == 0
    # invalid symbols are counted, read as N, never dropped
    p = parse_sequence("ACGT?Z")
    @test p.n_invalid == 2
    @test length(p.mask) == 6
    @test MitoASV.gc_content("GGCCAT") ≈ 4 / 6
    @test MitoASV.gc_content("NNNN") === missing
    @test MitoASV.gc_content("GSWA") ≈ 0.5
end

@testset "Stops, frames and partial sequences" begin
    orf = simulated_orf(60; seed = 3)
    @test !occursin('*', translate_sequence(orf))
    # internal stop
    bad = setcodon(orf, 20, "TAA")
    v = coding_validation(bad; frame = 1)
    @test v.n_internal_stops == 1
    @test v.record.internal_stop_positions == [20]
    # terminal stop is not internal; can be disallowed
    term = orf * "TAA"
    @test coding_validation(term; frame = 1).n_internal_stops == 0
    @test coding_validation(term; frame = 1).terminal_stop
    @test coding_validation(term; frame = 1, terminal_stop_allowed = false).n_internal_stops == 1
    # partial: incomplete terminal codon is not a frameshift
    part = orf * "AT"
    vp = coding_validation(part; frame = 1)
    @test vp.trailing_bases == 2
    @test vp.n_internal_stops == 0
    # leading partial codon with explicit frame
    lead = "GA" * orf
    @test translate_sequence(lead; frame = 3) == translate_sequence(orf)
    # unknown frame: frames evaluated, none chosen
    cv = coding_validation(orf)
    @test cv.frame === missing
    @test nrow(cv.frame_table) == 3
    @test cv.frame_table.open[1]
    ft = frame_candidates("ATGATGATG")
    @test ft.open == [true, true, true]
    @test MitoASV.frame_status(MitoASV.frame_stop_counts(MitoASV.parse_sequence("ATGATGATG").mask,
                                                          genetic_code(5))) == :multiple_open_frames
    # too short for a codon
    @test coding_validation("AT"; frame = 1).record.n_codons == 0
    # aligned input: full gap codon and partial gap codon
    al = coding_validation("ATG---TT-GGG"; frame = 1)
    @test al.protein == "M-XG"
    @test al.record.n_gap_codons == 1
    @test al.record.n_partial_gap_codons == 1
    @test al.length_nt == 8
end

@testset "GC by codon position and composition" begin
    r = coding_validation("GCAGCAGCA"; frame = 1)   # G C A ×3
    @test r.gc1 == 1.0 && r.gc2 == 1.0 && r.gc3 == 0.0
    comp = aa_composition("AAWX*-")
    @test comp[1] ≈ 2 / 3
    @test comp[19] ≈ 1 / 3       # W
    @test all(ismissing, aa_composition("X*-"))
end
