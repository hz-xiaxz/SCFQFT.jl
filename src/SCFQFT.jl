module SCFQFT

# Explicit imports from each package
using CompositeGrids: SimpleGrid
using GreenFunc.MeshGrids: FERMION
using GreenFunc.MeshArrays: MeshArray
using GreenFunc: MeshGrids
using SpecialFunctions: gamma

# Default mesh parameters
const DEFAULT_β = 0.1
const DEFAULT_Λ = 2.0
const DEFAULT_MESH_POINTS = 50

# Mesh creation function
function create_meshes(;
    β::Float64 = DEFAULT_β,
    Λ::Float64 = DEFAULT_Λ,
    n_points::Int = DEFAULT_MESH_POINTS,
)
    ω_mesh = MeshGrids.ImFreq(β, FERMION)
    k_mesh = SimpleGrid.Uniform([-sqrt(Λ), sqrt(Λ)], n_points)
    θ_mesh = SimpleGrid.Uniform([0, π], n_points)
    ϕ_mesh = SimpleGrid.Uniform([0, 2π], n_points)
    return ω_mesh, k_mesh, θ_mesh, ϕ_mesh
end

# Global meshes for default usage
const ω_mesh, k_mesh, θ_mesh, ϕ_mesh = create_meshes()

# Constants
const gamma_mat = [1 0; 0 -1]
const deltam = [1 0; 0 1]
# the cutoff constant, V should go to 0- for the renormalization procedure
const d = 3
const c = 2 * π^(d / 2) / (gamma(d / 2) * (2π)^d) * DEFAULT_Λ^(d - 2) / (d - 2)

# Helper functions for energy calculations
function ϵ_bar(; ϵ::Float64, v::Float64)
    ϵ + 4 / (3π * v)
end

function E(; ϵ::Float64, v::Float64, μ::Float64, Δ::Float64)
    √((ϵ_bar(ϵ = ϵ, v = v) - μ)^2 + Δ^2)
end

# Helper functions for quantum factors
function usq(; ϵ::Float64, v::Float64, μ::Float64, Δ::Float64)
    1 / 2 * (1 + (ϵ_bar(ϵ = ϵ, v = v) - μ) / (E(ϵ = ϵ, v = v, μ = μ, Δ = Δ) - μ))
end

function vsq(; ϵ::Float64, v::Float64, μ::Float64, Δ::Float64)
    1 / 2 * (1 - (ϵ_bar(ϵ = ϵ, v = v) - μ) / (E(ϵ = ϵ, v = v, μ = μ, Δ = Δ) - μ))
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
function G_mean(; ω::Float64, ϵ::Float64, v::Float64, μ::Float64, Δ::Float64)
    usq(ϵ = ϵ, v = v, μ = μ, Δ = Δ) / (-im * ω + E(ϵ = ϵ, v = v, μ = μ, Δ = Δ) - μ) +
    vsq(ϵ = ϵ, v = v, μ = μ, Δ = Δ) / (im * ω + E(ϵ = ϵ, v = v, μ = μ, Δ = Δ) - μ)
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
function F_mean(; ω::Float64, ϵ::Float64, v::Float64, μ::Float64, Δ::Float64)
    -Δ / (2 * (E(ϵ = ϵ, v = v, μ = μ, Δ = Δ) - μ)) * (
        1 / (-im * ω + E(ϵ = ϵ, v = v, μ = μ, Δ = Δ) - μ) +
        1 / (im * ω + E(ϵ = ϵ, v = v, μ = μ, Δ = Δ) - μ)
    )
end
# above ref Hausmann 2003 (3.63)-(3.67)

function Green(;
    v::Float64,
    μ::Float64,
    Δ::Float64,
    meshes = (ω_mesh, k_mesh, θ_mesh, ϕ_mesh),
)
    ω_m, k_m, θ_m, ϕ_m = meshes
    G_n = MeshArray(1:2, 1:2, k_m, θ_m, ϕ_m, ω_m; dtype = ComplexF64)
    for ind in eachindex(G_n)
        k = G_n.mesh[3][ind[3]]
        ω_n = G_n.mesh[6][ind[6]]
        ksq = k^2
        if ind[1] == 1 && ind[2] == 1
            G_n[ind] = G_mean(ω = ω_n, ϵ = ksq, v = v, μ = μ, Δ = Δ)
        elseif ind[1] == 1 && ind[2] == 2
            G_n[ind] = F_mean(ω = ω_n, ϵ = ksq, v = v, μ = μ, Δ = Δ)
        elseif ind[1] == 2 && ind[2] == 1
            G_n[ind] = -conj(F_mean(ω = -ω_n, ϵ = ksq, v = v, μ = μ, Δ = Δ))
        elseif ind[1] == 2 && ind[2] == 2
            G_n[ind] = -G_mean(ω = -ω_n, ϵ = ksq, v = v, μ = μ, Δ = Δ)
        end
    end
    return G_n
end
# compute M(K, Ω_n)
function SCF()
    # start from phase space component of G and F
    # the green function has no spin index, but Nambu indices
end

end # module
