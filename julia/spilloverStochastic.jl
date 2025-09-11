<< << << < HEAD
using ModelingToolkit, StochasticDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
using Plots

@parameters begin
    K = 1500, [description = "Maximum carrying capacity"]
    Kₛ = 0.85, [description = "Seasonal decrease of carrying capacity"]
    ϕ = 50, [description = "Shedding rate"]
    h = 0.05, [description = "Hunting efficiency"]
    hₛ = 0.1, [description = "Seasonal increase"]
    r = 0.25, [description = "Seasonal decrease"]
    μ = 0.05, [description = "Daily mortality"]
    ν = 0.7, [description = "Slower pathogen decay"]
    w_I = 0.25, [description = "Injury pathway weight"]
    w_P = 0.25, [description = "Processing pathway weight"]
    w_E = 0.25, [description = "Exposure pathway weight"]
    ρ_I = 0.05, [description = "Probability of transmission by injury"]
    ρ_P = 0.03, [description = "Probability of transmission by processing"]
    ρ_E = 0.05, [description = "Probability of transmission by eating"]
    l_I = 0.006, [description = "Probability injury occurrence"]
    l_P = 0.001, [description = "Probability processing occurrence"]
    l_E = 0.02, [description = "Probability eating occurrence"]
    η = 1e-4
    T = 4 / 24, [description = "Average daily hunting time"]
    Tₘ = 6 / 24, [description = "Increased hunting time"]
    λ = 0.08, [description = "Baseline transmission rate"]
    λᵦ = 1.2, [description = "Aggregation effect by precipitation"]
    k_c = ρ_P * w_P * l_P + ρ_E * w_E * l_E
    k_t = ρ_I * w_I * l_I
    # La Niña effect
    b = 1, [description = "Baseline precipitation"]
    ω = 3, [description = "Seasonal variation amplitude"]
    ϵ = 0.5, [description = "Noise intensity"]
    N = 2, [description = "La Niña effect intensity"]
    f = 5, [description = "Flood intensity"]
    Ω = 0.4, [description = "Optimal climate"]
end

@variables begin
    S_h(t)
    C_h(t)
    I_h(t)
    A(t)
    I_a(t)
    P(t)
    Γ_c(t)
    Γ_d(t)
    # Precipitation
    ξ(t)
    Λ(t)
    Kₑ(t)
    Tₑ(t)
    hₑ(t)
end

@brownians B

eqs = [
    ξ ~ abs((1 / (2 * ω)) * (b + ω * (1 + sin(2 * π * f * t))) - Ω) + B,
    Λ ~ λ * (1 + ξ * λᵦ),
    Kₑ ~ K * (1 - ξ * (1 - Kₛ)),
    Tₑ ~ (rem(t, 7) >= 5) * (((1 + 0.5 * ξ) / 24) * (Tₘ - T) + T) * ((ξ > 0.6) * ξ),
    hₑ ~ (rem(t, 7) >= 5) * h * (1 + hₛ * (ξ + Ω)),
    Γ_c ~ Tₑ * hₑ * k_c * (I_a / (I_a + A)),
    Γ_d ~ Tₑ * hₑ * k_t * (I_a / (I_a + A)),
    D(S_h) ~ -S_h * (1 - exp(Γ_d)),
    D(C_h) ~ -C_h * (1 - exp(Γ_c)),
    D(I_h) ~ S_h * (1 - exp(Γ_d)),
    D(A) ~ A * (r * (1 - (A + I_a) / Kₑ) - h - μ - Λ * (I_a / (I_a + A))),
    D(I_a) ~ Λ * A * (I_a / (I_a + A)) - I_a * (h + μ),
    D(P) ~ ϕ * I_a * (1 - I_a / Kₑ) - ν * P
]

@mtkcompile sp = System(eqs, t)
equations(sp)

u0 = [sp.S_h => 100.0,
    sp.I_h => 0.0,
    sp.A => 800.0,
    sp.I_a => 100.0,
    sp.P => 0.0,
    sp.C_h => 400.0]

prob = SDEProblem(sp, u0, (0.0, 600.0))
sol = solve(prob, SRIW1())

plot(sol; idxs=[I_h])
