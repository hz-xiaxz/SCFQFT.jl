using SCFQFT: ϵ_bar, E, usq, vsq, G_mean, F_mean
using Test

@testset "SCFQFT.jl" begin
    @testset "Energy calculations" begin
        # Test ϵ_bar
        @test SCFQFT.ϵ_bar(2.0) ≈ 1.0 + 4 / (3π * 2.0)

        # Test E with known values
        test_E = SCFQFT.E(v = 2.0, μ = 0.5, Δ = 0.1)
        @test test_E > 0  # Energy should be positive
        @test test_E isa Float64
    end

    @testset "Quantum factors" begin
        # Test usq + vsq = 1 (normalization)
        v, μ, Δ = 2.0, 0.5, 0.1
        @test SCFQFT.usq(v = v, μ = μ, Δ = Δ) + SCFQFT.vsq(v = v, μ = μ, Δ = Δ) ≈ 1.0
    end

    @testset "Green's functions" begin
        # Test G_mean and F_mean return expected types
        v, μ, Δ = 2.0, 0.5, 0.1
        @test G_mean(ω = 1.0, v = v, μ = μ, Δ = Δ) isa Complex{Float64}
        @test F_mean(ω = 1.0, v = v, μ = μ, Δ = Δ) isa Complex{Float64}

        # Test specific physical properties
        # Example: G_mean should be conjugate symmetric
        @test G_mean(ω = -1.0, v = v, μ = μ, Δ = Δ) ≈
              conj(G_mean(ω = 1.0, v = v, μ = μ, Δ = Δ))
    end
end
