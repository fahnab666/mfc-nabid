# CFD Equation Reference — Compressible Reactive Flow

Exact governing relations for MFC's domain. Use these verbatim — do not re-derive
from memory (mis-remembered signs and factors are silent physics bugs). Notation:
ρ density, u velocity, p pressure, E total energy per volume, e specific internal
energy, c sound speed, γ ratio of specific heats, Γ Grüneisen coefficient,
Y_k mass fraction, α_k volume fraction. Companion reasoning protocol:
`senior-cfd-physics.md`.

## 1. Compressible Euler equations (conservative form, what MFC solves)

    ∂ρ/∂t      + ∇·(ρu)          = 0
    ∂(ρu)/∂t   + ∇·(ρu⊗u + pI)   = 0
    ∂E/∂t      + ∇·((E + p)u)    = S_E     (reaction sources, P3 protocol)

with E = ρe + ½ρ|u|², closed by p = p(ρ, e, Y). Viscous terms (Navier-Stokes)
add ∇·τ and ∇·(τ·u − q) with τ = μ(∇u + ∇uᵀ) − (2/3)μ(∇·u)I — MFC's detonation
work is inviscid.

## 2. Five-equation multiphase model (Allaire/Kapila form)

    ∂(α_k ρ_k)/∂t + ∇·(α_k ρ_k u) = 0          (k = 1..N, one per fluid)
    ∂(ρu)/∂t      + ∇·(ρu⊗u + pI) = 0
    ∂E/∂t         + ∇·((E+p)u)    = 0
    ∂α_k/∂t       + u·∇α_k        = 0          (NON-conservative advection)

ρ = Σ α_k ρ_k, Y_k = α_k ρ_k / ρ. Single pressure and velocity; mixture-cell
closure is model-defined (Rocflu interpolation for JWL). Kapila's variant adds
K∇·u to the α equation — MFC's Allaire form does NOT.

## 3. Equations of state and sound speeds

Ideal gas:      p = (γ−1)ρe,             c² = γp/ρ,          Γ = γ−1
Stiffened gas:  p = (γ−1)ρe − γπ∞,       c² = γ(p+π∞)/ρ
  (Le Métayer convention; MFC namelist stores gamma = 1/(γ−1), pi_inf = γπ∞/(γ−1)… 
   see m_jwl.fpp init for the exact conversion — never guess the convention.)
Mie-Grüneisen:  p = p_ref(ρ) + Γρ(e − e_ref(ρ)),   Γ = V(∂p/∂e)|_V
JWL (products): p = A(1 − ω/(R₁V))e^(−R₁V) + B(1 − ω/(R₂V))e^(−R₂V) + ωρe,
                V = ρ₀/ρ
General frozen sound speed (any EOS):
    c² = (∂p/∂ρ)|_e + (p/ρ²)(∂p/∂e)|_ρ
Differentiate THROUGH every state-dependent coefficient (P1.3).

## 4. Rankine-Hugoniot jump conditions (frame of the shock, speed D)

    ρ₁(D − u₁) = ρ₂(D − u₂)                                (mass)
    p₂ − p₁ = ρ₁(D − u₁)(u₂ − u₁)                          (momentum)
    e₂ − e₁ = ½(p₁ + p₂)(v₁ − v₂),  v = 1/ρ               (energy; Hugoniot)
Rayleigh line: p₂ − p₁ = (ρ₁(D − u₁))²(v₁ − v₂)

Ideal-gas normal shock (Mach M₁ = (D−u₁)/c₁):
    p₂/p₁ = 1 + 2γ(M₁² − 1)/(γ+1)
    ρ₂/ρ₁ = (γ+1)M₁² / ((γ−1)M₁² + 2)      → (γ+1)/(γ−1) as M₁→∞
Strong-shock limit density ratio: 6 for γ=1.4.

## 5. Isentropic relations and Riemann invariants (ideal gas)

    p/ρ^γ = const;  T₂/T₁ = (p₂/p₁)^((γ−1)/γ)
    Riemann invariants: J± = u ± 2c/(γ−1) along dx/dt = u ± c
    Stagnation: p₀/p = (1 + (γ−1)/2 M²)^(γ/(γ−1))

## 6. Chapman-Jouguet and ZND detonation theory

CJ condition: the Rayleigh line is TANGENT to the fully-reacted Hugoniot; flow is
sonic in the frame of the front at the CJ plane: D − u_CJ = c_CJ.

With heat release q and effective γ of products (strong-detonation approx):
    D_CJ ≈ sqrt(2(γ² − 1)q)
    p_CJ ≈ ρ₀ D_CJ² / (γ + 1)
    ρ_CJ/ρ₀ ≈ (γ + 1)/γ,   u_CJ ≈ D_CJ/(γ + 1)
von Neumann spike (leading inert shock): p_VN ≈ 2 p_CJ (same-γ estimate).
ZND structure: inert shock to VN state → finite-rate reaction zone along the
Rayleigh line → sonic CJ plane → Taylor rarefaction (p drops to ≈0.35–0.40 p_CJ
at the center/piston). Reaction-zone length sets the resolution requirement.

## 7. Blast and explosion scaling

Sedov-Taylor strong point blast:  r_s(t) = ξ₀ (E t²/ρ₀)^(1/5)  (ξ₀ ≈ 1.03 air,
spherical); shock decays as p_s ∝ r_s⁻³ near-field.
Hopkinson-Cranz scaled distance: Z = R / W^(1/3)  (R in m, W in kg TNT) — all
far-field overpressure correlations (Kinney-Graham) are functions of Z only.
Acoustic far field: Δp ∝ 1/R.
Underwater (Cole): p_peak = K (W^(1/3)/R)^1.13, exponential decay time
θ ∝ W^(1/3)(W^(1/3)/R)^(−0.22).

## 8. Dimensionless groups (compute BEFORE choosing a model)

    Ma = u/c            compressibility (>0.3 ⇒ compressible)
    Re = ρuL/μ          viscous scale separation
    St = τ_p/τ_f        particle response (τ_p = ρ_p d²/(18μ) Stokes)
    Da = τ_f/τ_chem     reaction-flow coupling (Da >> 1 ⇒ thin front)
    CFL = (|u|+c)Δt/Δx  explicit stability (SSP-RK3 + upwind: CFL ≲ 1; MFC
                        default targets are tighter — check case defaults)

## 9. Particle equation of motion (dilute, point-particle)

    m_p du_p/dt = F_qs + F_pg + F_am
    F_qs = 3πμd C_D-corr (u_f − u_p)          quasi-steady drag (Re_p, Ma_p corr.)
    F_pg = V_p ρ_f Du_f/Dt                    pressure-gradient (dominant in blasts)
    F_am = ½ V_p ρ_f (Du_f/Dt − du_p/dt)      added mass (dominant for bubbles)
Energy: m_p c_p dT_p/dt = πd k_g Nu (T_f − T_p). Every gain here is a loss in the
gas-phase RHS (conservation pairing, P9/P3).

## 10. Discrete conservation identity (the property to protect)

A finite-volume update Q_i^{n+1} = Q_i^n − (Δt/Δx)(F_{i+1/2} − F_{i−1/2}) + Δt S_i
telescopes: Σ_i Q_i changes ONLY by boundary fluxes and Σ S_i. Any code path that
writes cell states directly (ghost fill, IBM correction, floors, clamps) sits
OUTSIDE this identity and is where conservation silently leaks — audit those
paths first when mass/energy drifts (see P8 Rankine-Hugoniot check).
