#==
based on https://docs.qojulia.org/examples/quantum-zeno-effect
==#
include("00 functions.jl")
using QuantumOptics, CairoMakie
# %%

# Parameters
Nc = 3
κ = 1.0
κ₂ = 0.3κ
g = 0.2κ
β = 0.2*sqrt(κ)
# gₛ = 5*1e-2κ
# T = [0:0.1:3/gₛ;]

# Basis and operators
bc = FockBasis(Nc)
ba = SpinBasis(1//2)

a = destroy(bc) ⊗ one(ba)
ad = create(bc) ⊗ one(ba)
σ⁺ = one(bc) ⊗ sigmap(ba)
σ = one(bc) ⊗ sigmam(ba)
n = σ⁺*σ

# Hamiltonian
H0(g) = g*ad*a*n
Hp(gₛ) = gₛ*(σ⁺ + σ)
Hf(β=β) = -1.0im*β*(ad - a) # Coherent drive of cavity

# Damping operators of master equation
J = [a]
rates = [κ]

# Stochastic damping operators
C = [sqrt(κ₂)*a]

# Initial state
ψ0 = fockstate(bc, 0) ⊗ spindown(ba)
ρ0 = dm(fockstate(bc, 0) ⊗ spindown(ba))

# Compute time evolution with no coupling to the cavity
free = timeevolution.master([0:0.1:π;], ρ0, Hp(1)+Hf(), J; rates=rates)


function zeno(;dt=1e-3, gₛ = 5*1e-2κ, T = [0:0.1:π/gₛ;])
	H = H0(g) + Hf(β) + Hp(gₛ)
	return stochastic.master(T, ρ0, H, J, C; rates=rates, dt=dt)
end

gₛ=[1e-1, 5e-2, 2.5e-2, 1e-2]
@time evolution = [zeno(gₛ=i) for i in gₛ];


proj = one(bc) ⊗ dm(spindown(ba)) # Calculate expectation value of atoms in |f>
p0(rho) = real(expect(proj, rho))

#%%
function plot_evolution(gs, ev)
	lines!(
		ev[1].*gs, p0(ev[2]),
		color=gs, colorrange=(0, maximum(gₛ)),
		label="$(gs)"
	)
end

f = Figure(size=fullsize)
ax = Axis(f[1,1], ylabel="<0|ρ|0>", xlabel="gₛt")
lines!(
	free[1], p0(free[2]),
	color=:black,
	label="1"
)
for (i, ev) in zip(gₛ, evolution)
	plot_evolution(i, ev)
end
axislegend(ax, position=:rb, framevisible=false, padding=0)
# Colorbar(f[2,1], colorrange=(0, maximum(gₛ)), vertical=false, label="gs", flipaxis=false)
f