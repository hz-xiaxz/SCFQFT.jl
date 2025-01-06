module SCFQFT
using LinearAlgebra
using GreenFunc, CompositeGrids

const beta = 0.1 # in the unit of T_F
const ω_mesh = MeshGrid.ImFreq(β, FERMION)
const kx_mesh = SimpleGrid.Uniform([-Λ, Λ], 50)
const ky_mesh = SimpleGrid.Uniform([-Λ, Λ], 50)
const kz_mesh = SimpleGrid.Uniform([-Λ, Λ], 50)

# Constants
const gamma = [1 0; 0 -1]
const deltam = [1 0; 0 1]
# the cutoff constant, V should go to 0- for the renormalization procedure
const Λ = 2
const d = 3
const c = 2 * π^(d / 2) / (SpecialFunctions.gamma(d / 2) * (2π)^d) * Λ^(d - 2) / (d - 2)

# Helper functions for energy calculations
function ϵ_bar(v::Float64)
    1 + 4 / (3π * v)
end

function E(; v::Float64, μ::Float64, Δ::Float64)
    √((ϵ_bar(v) - μ)^2 + Δ^2)
end

# Helper functions for quantum factors
function usq(; v::Float64, μ::Float64, Δ::Float64)
    1 / 2 * (1 + (ϵ_bar(v) - μ) / (E(v = v, μ = μ, Δ = Δ) - μ))
end

function vsq(; v::Float64, μ::Float64, Δ::Float64)
    1 / 2 * (1 - (ϵ_bar(v) - μ) / (E(v = v, μ = μ, Δ = Δ) - μ))
end

# Main exported functions
"""
Calculate the mean Green's function G
Parameters:
    ω::Float64 - frequency
    ϵ::Float64 - energy
    v::Float64 - coupling strength
    μ::Float64 - chemical potential
    Δ::Float64 - gap parameter
"""
function G_mean(; ω::Float64, v::Float64, μ::Float64, Δ::Float64)
    usq(v = v, μ = μ, Δ = Δ) / (-im * ω + E(v = v, μ = μ, Δ = Δ) - μ) +
    vsq(v = v, μ = μ, Δ = Δ) / (im * ω + E(v = v, μ = μ, Δ = Δ) - μ)
end

"""
Calculate the mean anomalous Green's function F
Parameters:
    ω::Float64 - frequency
    ϵ::Float64 - energy
    v::Float64 - coupling strength
    μ::Float64 - chemical potential
    Δ::Float64 - gap parameter
"""
function F_mean(; ω::Float64, v::Float64, μ::Float64, Δ::Float64)
    -Δ / (2 * (E(v = v, μ = μ, Δ = Δ) - μ)) *
    (1 / (-im * ω + E(v = v, μ = μ, Δ = Δ) - μ) + 1 / (im * ω + E(v = v, μ = μ, Δ = Δ) - μ))
end
# above ref Hausmann 2003 (3.63)-(3.67)

function Green(; v::Float64, μ::Float64, Δ::Float64)
    G_n = MeshArray(1:2, 1:2, kx_mesh, ky_mesh, kz_mesh, ω_mesh; dtype = ComplexF64)
    for ind in eachindex(G_n)
        kx = G_n.mesh[3][ind[3]]
        ky = G_n.mesh[4][ind[4]]
        kz = G_n.mesh[5][ind[5]]
        ω_n = G_n.mesh[6][ind[6]]
        if ind[1] == 1 && ind[2] == 1
            G_n[ind] = G_mean(ω = ω_n, v = v, μ = μ, Δ = Δ)
        elseif ind[1] == 1 && ind[2] == 2

        end
    end
end
# compute M(K, Ω_n)
function SCF()
    # start from phase space component of G and F
    # the green function has no spin index, but Nambu indices
end

end # module
