using StochasticDiffEq
using LinearAlgebra
using Plots

# Parameters

K = 1500
Kₛ = 0.85
ϕ = 50
h = 0.05
hₛ = 0.1
r = 0.25
μ = 0.05
ν = 0.7
w_I = 0.25
w_P = 0.25
w_E = 0.25
ρ_I = 0.05
ρ_P = 0.03
ρ_E = 0.05
l_I = 0.006
l_P = 0.001
l_E = 0.02
η = 1e-4
T = 4 / 24
Tₘ = 6 / 24
λ = 0.08
λᵦ = 1.2
k_c = ρ_P * w_P * l_P + ρ_E * w_E * l_E
k_t = ρ_I * w_I * l_I
# La Niña effect
b = 1
ω = 3
ϵ = 100
N = 2
f = 5
Ω = 0

function system!(du, u, p, t)
    du[1] = abs((1 / (2 * ω)) * (b + ω * (1 + sin(2 * π * (t / 365) - π / 2))) - Ω)
    du[2] = λ * (1 + du[1] * λᵦ)
    du[3] = K * (1 - du[1] * (1 - Kₛ))
    du[4] = (rem(t, 7) >= 5) * (((1 + 0.5 * du[1]) / 24) * (Tₘ - T) + T) * ((du[1] > 0.6) * du[1])
    du[5] = (rem(t, 7) >= 5) * h * (1 + hₛ * (du[1] + Ω))
    du[6] = du[4] * du[5] * k_c * (u[4] / (u[4] + u[3]))
    du[7] = du[4] * du[5] * k_t * (u[4] / (u[4] + u[3]))
    du[8] = -u[1] * (1 - exp(du[7]))
    du[9] = -u[6] * (1 - exp(du[6]))
    du[10] = u[1] * (1 - exp(du[7]))
    du[11] = u[3] * (r * (1 - (u[3] + u[4]) / du[3]) - h - μ - du[2] * (u[4] / (u[4] + u[3])))
    du[12] = du[2] * u[3] * (u[4] / (u[4] + u[3])) - u[4] * (h + μ)
    du[13] = ϕ * u[4] * (1 - u[4] / du[3]) - ν * u[5]
end

function noise!(du, u, p, t)
    @. du = 0.1
end

M = Diagonal(ones(13))

u0 = [100.0, 0.0, 800.0, 100.0, 0.0, 400.0, zeros(7)...]
prob = SDEProblem(SDEFunction(system!, noise!; mass_matrix=M), u0, (0.0, 10.0))
sol = solve(prob, ImplicitEM())

plot(sol)
