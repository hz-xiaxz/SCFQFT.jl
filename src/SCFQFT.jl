module SCFQFT

# Explicit imports from each package
using CompositeGrids: SimpleGrid
using GreenFunc.MeshGrids: FERMION, BOSON
using GreenFunc.MeshArrays: MeshArray
using GreenFunc: MeshGrids
using SpecialFunctions: gamma
using LinearAlgebra: norm

# Default mesh parameters
const DEFAULT_β = 0.1
const DEFAULT_Λ = 2.0
const DEFAULT_MESH_POINTS = 50

# Define a new struct to hold the variables
@kwdef struct Parameters
    v::Float64
    μ::Float64
    Δ::Float64
end

# Mesh creation function
function create_meshes(;
    β::Float64 = DEFAULT_β,
    Λ::Float64 = DEFAULT_Λ,
    n_points::Int = DEFAULT_MESH_POINTS,
)
    ω_mesh = MeshGrids.ImFreq(β, FERMION)
    Ω_mesh = MeshGrids.ImFreq(β, BOSON)
    k_mesh = SimpleGrid.Uniform([-sqrt(Λ), sqrt(Λ)], n_points)
    θ_mesh = SimpleGrid.Uniform([0, π], n_points)
    ϕ_mesh = SimpleGrid.Uniform([0, 2π], n_points)
    return ω_mesh, Ω_mesh, k_mesh, θ_mesh, ϕ_mesh
end

# Global meshes for default usage
const ω_mesh, Ω_mesh, k_mesh, θ_mesh, ϕ_mesh = create_meshes()

# Constants
const γ = [1 0; 0 -1]

# the cutoff constant, V should go to 0- for the renormalization procedure
const d = 3
const c = 2 * π^(d / 2) / (gamma(d / 2) * (2π)^d) * DEFAULT_Λ^(d - 2) / (d - 2)

# Helper functions for energy calculations
function ϵ_bar(; ϵ::Float64, v::Float64)
    ϵ + 4 / (3π * v)
end

function E(; ϵ::Float64, para::Parameters)
    √((ϵ_bar(ϵ = ϵ, v = para.v) - para.μ)^2 + para.Δ^2)
end

# Helper functions for quantum factors
function usq(; ϵ::Float64, para::Parameters)
    1 / 2 * (1 + (ϵ_bar(ϵ = ϵ, v = para.v) - para.μ) / (E(ϵ = ϵ, para = para) - para.μ))
end

function vsq(; ϵ::Float64, para::Parameters)
    1 / 2 * (1 - (ϵ_bar(ϵ = ϵ, v = para.v) - para.μ) / (E(ϵ = ϵ, para = para) - para.μ))
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
function G_mean(; ω::Float64, ϵ::Float64, para::Parameters)
    usq(ϵ = ϵ, para = para) / (-im * ω + E(ϵ = ϵ, para = para) - para.μ) -
    vsq(ϵ = ϵ, para = para) / (im * ω + E(ϵ = ϵ, para = para) - para.μ)
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
function F_mean(; ω::Float64, ϵ::Float64, para::Parameters)
    -para.Δ / (2 * (E(ϵ = ϵ, para = para) - para.μ)) * (
        1 / (-im * ω + E(ϵ = ϵ, para = para) - para.μ) +
        1 / (im * ω + E(ϵ = ϵ, para = para) - para.μ)
    )
end
# above ref Hausmann 2003 (3.63)-(3.67)

function G_n(; para::Parameters, meshes = (ω_mesh, k_mesh, θ_mesh, ϕ_mesh))
    ω_m, k_m, θ_m, ϕ_m = meshes
    G_n = MeshArray(1:2, 1:2, k_m, θ_m, ϕ_m, ω_m; dtype = ComplexF64)
    for ind in eachindex(G_n)
        k = G_n.mesh[3][ind[3]]
        ω_n = G_n.mesh[6][ind[6]]
        ksq = k^2
        if ind[1] == 1 && ind[2] == 1
            G_n[ind] = G_mean(ω = ω_n, ϵ = ksq, para = para)
        elseif ind[1] == 1 && ind[2] == 2
            G_n[ind] = F_mean(ω = ω_n, ϵ = ksq, para = para)
        elseif ind[1] == 2 && ind[2] == 1
            G_n[ind] = -conj(F_mean(ω = -ω_n, ϵ = ksq, para = para))
        elseif ind[1] == 2 && ind[2] == 2
            G_n[ind] = -G_mean(ω = -ω_n, ϵ = ksq, para = para)
        end
    end
    return G_n
end

"""Calculate k·K term in spherical coordinates"""
function calc_dot_product(
    K::Float64,
    k::Float64,
    Θ::Float64,
    θ::Float64,
    Φ::Float64,
    ϕ::Float64,
)
    K * k * (sin(Θ) * sin(θ) * cos(Φ - ϕ) + cos(Θ) * cos(θ))
end

"""Calculate k±K/2 squared terms"""
function calc_k_squared_terms(K::Float64, k::Float64, dot_term::Float64)
    square_term = K^2 / 4 + k^2
    ksq1 = square_term - dot_term
    ksq2 = square_term + dot_term
    return ksq1, ksq2
end

"""Calculate frequency sum for given Nambu indices"""
function calc_frequency_sum(
    α1::Int,
    α2::Int,
    Ω_n::Float64,
    ksq1::Float64,
    ksq2::Float64,
    ω_m::AbstractArray,
    para::Parameters,
)

    omega_sum = 0.0im
    @inbounds for ω_n in ω_m
        if α1 == 1 && α2 == 1
            omega_sum +=
                G_mean(ω = Ω_n - ω_n, ϵ = ksq1, para = para) *
                G_mean(ω = ω_n, ϵ = ksq2, para = para)
        elseif α1 == 1 && α2 == 2
            omega_sum +=
                F_mean(ω = Ω_n - ω_n, ϵ = ksq1, para = para) *
                F_mean(ω = ω_n, ϵ = ksq2, para = para)
        elseif α1 == 2 && α2 == 1
            omega_sum +=
                conj(F_mean(ω = -(Ω_n - ω_n), ϵ = ksq1, para = para)) *
                conj(F_mean(ω = -ω_n, ϵ = ksq2, para = para))
        else # α1 == 2 && α2 == 2
            omega_sum +=
                G_mean(ω = -(Ω_n - ω_n), ϵ = ksq1, para = para) *
                G_mean(ω = -ω_n, ϵ = ksq2, para = para)
        end
    end
    return omega_sum
end

function M_n_atomic(;
    para::Parameters,
    meshes = (ω_mesh, Ω_mesh, k_mesh, θ_mesh, ϕ_mesh),
    ΔV::Float64,
    K::Float64,
    Θ::Float64,
    Φ::Float64,
    Ω_n::Float64,
    α1::Int,
    α2::Int,
)
    ω_m, Ω_m, k_m, θ_m, ϕ_m = meshes

    k_sum = 0.0im
    @inbounds for k_idx in eachindex(k_m), θ_idx in eachindex(θ_m), ϕ_idx in eachindex(ϕ_m)

        k = k_m[k_idx]
        θ = θ_m[θ_idx]
        ϕ = ϕ_m[ϕ_idx]

        dot_term = calc_dot_product(K, k, Θ, θ, Φ, ϕ)
        ksq1, ksq2 = calc_k_squared_terms(K, k, dot_term)

        omega_sum = calc_frequency_sum(α1, α2, Ω_n, ksq1, ksq2, ω_m, para)
        k_sum += omega_sum / DEFAULT_β * k^2 * sin(θ) / (2π)^3 * ΔV

        if α1 == α2
            k_sum -= 1 / 2 * ΔV / (2π)^3 * sin(Θ)
            # 1/2 is due to the unit of E_F
        end
    end
    return k_sum
end

"""
The non-divergent part of the pair propagator ``\\chi``
"""
function M_n(; para::Parameters, meshes = (ω_mesh, Ω_mesh, k_mesh, θ_mesh, ϕ_mesh))
    ω_m, Ω_m, k_m, θ_m, ϕ_m = meshes
    M_n = MeshArray(1:2, 1:2, k_m, θ_m, ϕ_m, Ω_m; dtype = ComplexF64)
    Δk = k_m[2] - k_m[1]
    Δθ = θ_m[2] - θ_m[1]
    Δϕ = ϕ_m[2] - ϕ_m[1]
    ΔV = Δk * Δθ * Δϕ

    @inbounds for ind in eachindex(M_n)
        K = k_m[ind[3]]
        Θ = θ_m[ind[4]]
        Φ = ϕ_m[ind[5]]
        Ω_n = Ω_m[ind[6]]

        M_n[ind] = M_n_atomic(
            para = para,
            meshes = meshes,
            ΔV = ΔV,
            K = K,
            Θ = Θ,
            Φ = Φ,
            Ω_n = Ω_n,
            α1 = ind[1],
            α2 = ind[2],
        )
    end
    return M_n
end

function Γ_n_atomic(; T::Float64, m_n::ComplexF64, α1::Int, α2::Int)
    if α1 == α2
        return inv(1 / T + m_n)
    else
        return 1 / m_n
    end
end

function Γ_n(; para::Parameters, meshes = (ω_mesh, Ω_mesh, k_mesh, θ_mesh, ϕ_mesh))
    ω_m, Ω_m, k_m, θ_m, ϕ_m = meshes
    m_n = M_n(para = para, meshes = meshes)
    T = para.v / (8π)
    Γ_n = MeshArray(1:2, 1:2, k_m, θ_m, ϕ_m, Ω_m; dtype = ComplexF64)
    @inbounds for ind in eachindex(Γ_n)
        Γ_n[ind] = Γ_n_atomic(T = T, m_n = m_n[ind], α1 = ind[1], α2 = ind[2])
    end
    return Γ_n
end

function Self_energy_atomic(;
    para::Parameters,
    meshes = (ω_mesh, Ω_mesh, k_mesh, θ_mesh, ϕ_mesh),
    α1::Int,
    α2::Int,
    k::Float64,
    θ::Float64,
    ϕ::Float64,
    ω::Float64,
)
    ω_m, Ω_m, k_m, θ_m, ϕ_m = meshes
    Δk = k_m[2] - k_m[1]
    Δθ = θ_m[2] - θ_m[1]
    Δϕ = ϕ_m[2] - ϕ_m[1]
    ΔV = Δk * Δθ * Δϕ
    Σ_sum = 0.0im
    if α1 == 1 && α2 == 2
        Σ_sum += para.Δ
    elseif α1 == 2 && α1 == 1
        Σ_sum += conj(para.Δ)
    end
    # this choice of k_1 meshes may cause some problems
    @inbounds for K in k_m, Θ in θ_m, Φ in ϕ_m
        ksqG = K^2
        K_vec = [K * sin(Θ) * cos(Φ), K * sin(Θ) * sin(Φ), K * cos(Θ)]
        k_vec = [k * sin(θ) * cos(ϕ), k * sin(θ) * sin(ϕ), k * cos(θ)]
        diff = K_vec - k_vec
        ω_sum = zero(Complex64)
        @inbounds for ω_n in ω_m
            G = zero(Complex64)
            if α2 == 1 && α1 == 1
                G = G_mean(ω = ω_n, ϵ = ksqG, para = para)
            elseif α2 == 1 && α1 == 2
                G = F_mean(ω = ω_n, ϵ = ksqG, para = para)
            elseif α2 == 2 && α1 == 1
                G = -conj(F_mean(ω = -ω_n, ϵ = ksqG, para = para))
            elseif α2 == 2 && α1 == 2
                G = -G_mean(ω = -ω_n, ϵ = ksqG, para = para)
            end
            diff_module = norm(diff)
            diff_θ = acos(diff[3] / diff_module)
            diff_ϕ = mod2pi(atan(diff[2], diff[1]))
            m = M_n_atomic(
                para = para,
                meshes = meshes,
                ΔV = ΔV,
                K = diff_module,
                Θ = diff_θ,
                Φ = diff_ϕ,
                Ω_n = ω - ω_n,
                α1 = α1,
                α2 = α2,
            )
            Γ = Γ_n_atomic(T = para.v / (8π), m_n = m, α1 = α1, α2 = α2)
            ω_sum += G * Γ
        end
        Σ_sum += ω_sum / DEFAULT_β * ΔV / (2π)^3 * K^2 * sin(Θ)
    end
    return Σ_sum
end


function Σ_n(; para::Parameters, meshes = (ω_mesh, Ω_mesh, k_mesh, θ_mesh, ϕ_mesh))
    ω_m, Ω_m, k_m, θ_m, ϕ_m = meshes
    Σ_n = MeshArray(1:2, 1:2, k_m, θ_m, ϕ_m, ω_m; dtype = ComplexF64)
    @inbounds for ind in eachindex(Σ_n)
        Σ_n[ind] = Self_energy_atomic(
            para = para,
            meshes = meshes,
            α1 = ind[1],
            α2 = ind[2],
            k = k_m[ind[3]],
            θ = θ_m[ind[4]],
            ϕ = ϕ_m[ind[5]],
            ω = ω_m[ind[6]],
        )
    end
    return Σ_n
end

function SCF()
    # start from phase space component of G and F
    # the green function has no spin index, but Nambu indices
end
end
