using Test
using MeanFieldHom
using TensND

@testset "Core — traits" begin
    C_iso = TensISO{3}(3.0, 2.0)
    @test MeanFieldHom.material_symmetry(C_iso) isa MeanFieldHom.IsotropicSym

    # Analytical / Residue / DECUHR singletons
    @test MeanFieldHom.Analytical() isa MeanFieldHom.AbstractAlgorithm
    @test MeanFieldHom.Residue() isa MeanFieldHom.AbstractAlgorithm
    @test MeanFieldHom.DECUHR() isa MeanFieldHom.AbstractAlgorithm
end
