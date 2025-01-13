using Test
using SCFQFT:
    calc_dot_product,
    calc_k_squared_terms,
    calc_frequency_sum,
    M_n,
    create_meshes,
    Parameters,
    DEFAULT_β,
    Γ_n
using GreenFunc: MeshGrids, FERMION, BOSON

@testset "M_n calculations" begin
    @testset "Dot product calculation" begin
        # Test orthogonal vectors
        @test calc_dot_product(1.0, 1.0, 0.0, Float64(π / 2), 0.0, 0.0) ≈ 0.0 atol = 1e-10

        # Test parallel vectors
        @test calc_dot_product(1.0, 1.0, 0.0, 0.0, 0.0, 0.0) ≈ 1.0 atol = 1e-10

        # Test antiparallel vectors
        @test calc_dot_product(1.0, 1.0, 0.0, Float64(π), 0.0, 0.0) ≈ -1.0 atol = 1e-10
    end

    @testset "k squared terms" begin
        # Test zero momentum
        ksq1, ksq2 = calc_k_squared_terms(0.0, 1.0, 0.0)
        @test ksq1 ≈ ksq2 ≈ 1.0

        # Test with finite K and dot product
        ksq1, ksq2 = calc_k_squared_terms(2.0, 1.0, 1.0)
        @test ksq1 + ksq2 ≈ 2 * (1.0 + 1.0) # sum should be independent of dot product
    end

    @testset "M_n calculations" begin
        # Create Parameters struct for testing
        params = Parameters(v = 1.0, μ = 0.5, Δ = 0.1)  # Use kwarg syntax for Parameters

        # Test M_n function with the parameters struct
        ω_mesh, Ω_mesh, k_mesh, θ_mesh, ϕ_mesh =
            create_meshes(β = 0.1, Λ = 2.0, n_points = 5)
        Ω_mesh = MeshGrids.ImFreq(0.1, BOSON; Euv = 1.0)
        ω_mesh = MeshGrids.ImFreq(0.1, FERMION; Euv = 1.0)
        M_result = M_n(para = params, meshes = (ω_mesh, Ω_mesh, k_mesh, θ_mesh, ϕ_mesh))

        @test !isnothing(M_result)  # Ensure M_n returns a result
        @test size(M_result) == (2, 2, 5, 5, 5, 36)  # Check dimensions with n_points=5
        @test eltype(M_result) == ComplexF64  # Check type
    end
end

@testset "Bethe-Salpeter equation" begin
    @testset "Γ_n Tests" begin
        # Setup test meshes
        ω_mesh = range(-1, 1, length = 3)
        Ω_mesh = range(-1, 1, length = 3)
        k_mesh = range(0, 1, length = 2)
        θ_mesh = range(0, π, length = 2)
        ϕ_mesh = range(0, 2π, length = 2)

        # Setup parameters
        para = Parameters(v = 8π, μ = 0.5, Δ = 0.1) # This should make T=1 for simpler testing

        # Calculate Γ_n
        result = Γ_n(para = para, meshes = (ω_mesh, Ω_mesh, k_mesh, θ_mesh, ϕ_mesh))

        # Test shape and type
        @test eltype(result) == ComplexF64

    end
end
