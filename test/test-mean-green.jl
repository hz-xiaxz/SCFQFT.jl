using SCFQFT: ϵ_bar, E, usq, vsq, G_mean, F_mean, G_n, create_meshes
using Test
using GreenFunc
using CompositeGrids

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
    @test SCFQFT.usq(ϵ = ϵ, v = v, μ = μ, Δ = Δ) + SCFQFT.vsq(ϵ = ϵ, v = v, μ = μ, Δ = Δ) ≈
          1.0
end

@testset "Green function" begin
    # Test parameters
    v = 1.0
    μ = 0.5
    Δ = 0.1

    # Create small test meshes
    β = 0.1
    Λ = 2.0
    n_points = 5  # Small number for testing
    test_meshes = create_meshes(β = β, Λ = Λ, n_points = n_points)
    ω_m, Ω_m, k_m, θ_m, ϕ_m = test_meshes

    # Get Green's function with test meshes
    G = G_n(v = v, μ = μ, Δ = Δ, meshes = (ω_m, k_m, θ_m, ϕ_m))

    @testset "Basic properties" begin
        # Test return value
        @test !isnothing(G)

        # Test dimensions (2×2 Nambu space × r × θ × φ × ω)
        expected_size = (2, 2, n_points, n_points, n_points, length(test_meshes[1]))
        @test size(G) == expected_size

        # Test type
        @test eltype(G) == ComplexF64
    end

    @testset "UV cutoff" begin
        # Test that Green function vanishes for r > Λ
        for r_idx = 1:n_points
            r = test_meshes[2][r_idx]  # k_mesh contains radial values
            if r > Λ
                for θ_idx = 1:n_points, ϕ_idx = 1:n_points
                    @test all(abs.(G[:, 1, r_idx, θ_idx, ϕ_idx, :]) .< 1e-10)
                    @test all(abs.(G[:, 2, r_idx, θ_idx, ϕ_idx, :]) .< 1e-10)
                end
            end
        end
    end

    @testset "Spherical symmetry" begin
        # Test that G only depends on |k| = r, not on angles
        r_idx, ω_idx = 1, 1
        reference_value = G[1, 1, r_idx, 1, 1, ω_idx]

        for θ_idx = 1:n_points, ϕ_idx = 1:n_points
            @test G[1, 1, r_idx, θ_idx, ϕ_idx, ω_idx] ≈ reference_value rtol = 1e-10
        end
    end
end
