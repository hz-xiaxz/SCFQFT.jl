using Test
using SCFQFT: Parameters, create_meshes, Self_energy_atomic, Σ_n
using LinearAlgebra
using GreenFunc
using CompositeGrids

@testset "SCFQFT.jl" begin
    # Setup default parameters and meshes
    para = Parameters(v = 1.0, μ = 0.5, Δ = 0.1)
    meshes = create_meshes(
        β = 0.1,
        Λ = 1.0,
        n_points = 2,  # Reduced for testing
    )
    ω_m, Ω_m, k_m, θ_m, ϕ_m = meshes
    ΔV = (k_m[2] - k_m[1]) * (θ_m[2] - θ_m[1]) * (ϕ_m[2] - ϕ_m[1])

    @testset "Self_energy_atomic" begin
        # Test diagonal elements
        Σ11 = Self_energy_atomic(
            para = para,
            meshes = meshes,
            α1 = 1,
            α2 = 1,
            k = 1.0,
            θ = π / 2,
            ϕ = 0.0,
            ω = 0.0,
        )
        @test isa(Σ11, Complex)
        @test !isnan(real(Σ11))
        @test !isnan(imag(Σ11))

        # Test anomalous components
        Σ12 = Self_energy_atomic(
            para = para,
            meshes = meshes,
            α1 = 1,
            α2 = 2,
            k = 1.0,
            θ = π / 2,
            ϕ = 0.0,
            ω = 0.0,
        )

        # Test Hermiticity
        Σ21 = Self_energy_atomic(
            para = para,
            meshes = meshes,
            α1 = 2,
            α2 = 1,
            k = 1.0,
            θ = π / 2,
            ϕ = 0.0,
            ω = 0.0,
        )
        @test isapprox(Σ21, conj(Σ12), rtol = 1e-10)
    end

    @testset "Full Self_energy" begin
        Σ = Σ_n(para = para, meshes = meshes)

        # Test dimensions
        expected_size = (2, 2, length(k_m), length(θ_m), length(ϕ_m), length(ω_m))
        @test size(Σ) == expected_size
        @test eltype(Σ) == ComplexF64

        # Test symmetries
        @inbounds for ind in eachindex(Σ)
            α1, α2, k_idx, θ_idx, ϕ_idx, ω_idx = Tuple(ind)

            # Hermiticity
            @test isapprox(
                Σ[2, 1, k_idx, θ_idx, ϕ_idx, ω_idx],
                conj(Σ[1, 2, k_idx, θ_idx, ϕ_idx, ω_idx]),
                rtol = 1e-10,
            )

            # Reality of trace
            @test abs(
                imag(
                    Σ[1, 1, k_idx, θ_idx, ϕ_idx, ω_idx] +
                    Σ[2, 2, k_idx, θ_idx, ϕ_idx, ω_idx],
                ),
            ) < 1e-10
        end

        # Test rotational invariance
        @inbounds for k_idx = 1:(length(k_m)÷2),
            θ_idx = 1:length(θ_m),
            ω_idx = 1:length(ω_m)

            @test isapprox(
                Σ[1, 1, k_idx, θ_idx, 1, ω_idx],
                Σ[1, 1, k_idx, θ_idx, end, ω_idx],
                rtol = 1e-10,
            )
        end

        # Test momentum scaling
        # k_mid = length(k_m) ÷ 2
        # low_k_norm = norm(view(Σ, :, :, 1, :, :, :))
        # mid_k_norm = norm(view(Σ, :, :, k_mid, :, :, :))
        # @test low_k_norm < mid_k_norm
    end
end
