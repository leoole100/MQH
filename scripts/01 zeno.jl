#==
based on https://docs.qojulia.org/examples/quantum-zeno-effect
==#
include("00 functions.jl")

# %%

# Parameters
Nc = 3
κ = 1.0
κ₂ = 0.3κ
g = 0.2κ
β = 0.2*sqrt(κ)
gₛ = 5*1e-2κ
T = [0:0.1:3/gₛ;]

# Basis and operators
bc = FockBasis(Nc)
ba = SpinBasis(1//2)

a = destroy(bc) ⊗ one(ba)
ad = create(bc) ⊗ one(ba)
σ⁺ = one(bc) ⊗ sigmap(ba)
σ = one(bc) ⊗ sigmam(ba)
n = σ⁺*σ

# %%

# Hamiltonian
H0 = g*ad*a*n
Hp = gₛ*(σ⁺ + σ)

# Coherent drive of cavity
Hf = -1.0im*β*(ad - a)

H = H0 + Hf + Hp

# Damping operators of master equation
J = [a]
rates = [κ]

# %%

# Stochastic damping operators
C = [sqrt(κ₂)*a]

# Initial state
ψ0 = fockstate(bc, 0) ⊗ spindown(ba)
ρ0 = dm(fockstate(bc, 0) ⊗ spindown(ba))

# Solve stochastic master equation
dt = 1e-3
@time tout, ρt = stochastic.master(T, ρ0, H, J, C; rates=rates, dt=dt);

# %%

C_zeno = [sqrt(κ)*a]
gₛ_zeno = 1e-3κ
Hp_zeno = gₛ_zeno*(σ⁺ + σ)
H_zeno = H0 + Hf + Hp_zeno
T_zeno = [0:1:3/gₛ_zeno;]

@time tout_zeno, ρt_zeno = stochastic.master(T_zeno, ρ0, H_zeno, J, C_zeno; rates=rates, dt=dt);

# %%
# Compute time evolution with no coupling to the cavity
@time tout2, ρt_det = timeevolution.master(T, ρ0, H - g*ad*a*n, J; rates=rates);

#%%
# Calculate expectation value of atoms in |f>
proj = one(bc) ⊗ dm(spindown(ba))
p0 = real(expect(proj, ρt))
p0_zeno = real(expect(proj, ρt_zeno))
p0_det = real(expect(proj, ρt_det))

# %%
f = Figure(size=fullsize)
a = Axis(f[1,1], ylabel="<0|ρ|0>", xlabel="gₛt")
lines!(tout.*gₛ, p0, label="low zeno")
lines!(tout_zeno.*gₛ_zeno, p0_zeno, label="High zeno")
lines!(tout2.*gₛ, p0_det, label="Free")
f