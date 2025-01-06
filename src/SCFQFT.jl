module SCFQFT
using LinearAlgebra


# Constants
const gamma = [1 0; 0 -1]
const deltam = [1 0; 0 1]

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

end # module
