using Test
using SCFQFT: calc_dot_product, calc_k_squared_terms, calc_frequency_sum, M_n

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

end
