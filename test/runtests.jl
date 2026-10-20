using MitoASV
using DataFrames
using Test
using Statistics

include("helpers.jl")

@testset "MitoASV" begin
    include("test_translation.jl")
    include("test_hydrophobicity.jl")
    include("test_interface.jl")
    include("test_substitutions.jl")
    include("test_codons.jl")
    include("test_clustering.jl")
    include("test_conservation.jl")
    include("test_groups.jl")
    include("test_anomaly.jl")
    include("test_concordance.jl")
    include("test_edgecases.jl")
    include("test_species.jl")
    include("test_genecolumns.jl")
end
