using SCFQFT: ϵ_bar, E, usq, vsq, G_mean, F_mean
using Test

@testset "SCFQFT.jl" begin
    @testset "Energy calculations" begin
        # Test ϵ_bar
        @test SCFQFT.ϵ_bar(ϵ = 1.0, v = 2.0) ≈ 1.0 + 4 / (3π * 2.0)

        # Test E with known values
        test_E = SCFQFT.E(ϵ = 1.0, v = 2.0, μ = 0.5, Δ = 0.1)
        @test test_E > 0  # Energy should be positive
        @test test_E isa Float64
    end

    @testset "Quantum factors" begin
        # Test usq + vsq = 1 (normalization)
        ϵ, v, μ, Δ = 1.0, 2.0, 0.5, 0.1
        @test SCFQFT.usq(ϵ = ϵ, v = v, μ = μ, Δ = Δ) +
              SCFQFT.vsq(ϵ = ϵ, v = v, μ = μ, Δ = Δ) ≈ 1.0
    end

    @testset "Green's functions" begin
        # Test G_mean and F_mean return expected types
        params = (ω = 1.0, ϵ = 1.0, v = 2.0, μ = 0.5, Δ = 0.1)
        @test G_mean(; params...) isa Complex{Float64}
        @test F_mean(; params...) isa Complex{Float64}

        # Test specific physical properties
        # Example: G_mean should be conjugate symmetric
        @test G_mean(ω = -1.0, ϵ = 1.0, v = 2.0, μ = 0.5, Δ = 0.1) ≈
              conj(G_mean(ω = 1.0, ϵ = 1.0, v = 2.0, μ = 0.5, Δ = 0.1))
    end
end
