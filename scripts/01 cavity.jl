#=
Implementation of a optomechanical cavity.
based on: https://docs.qojulia.org/examples/optomech-cooling/

$$
	H = -\Delta a^\dagger a + \omega_m b^\dagger b + g (a^\dagger a)(b^\dagger + b) + \eta (a + a^\dagger)
$$

where $a$ and $b$ are the annihilation operators for the cavity and mechanical mode, respectively. The cavity is driven by a coherent drive with amplitude $\eta$ and detuning $\Delta$. The mechanical mode has frequency $\omega_m$ and is coupled to the cavity with strength $g$.
=#

using QuantumOptics
using CairoMakie

#%%

# Parameters
ω_mech = 1.
Δ = -ω_mech

# Constants
g = 1.
η = 1.
κ = 0.1	# cavity decay rate
γ = 0.1	# mechanical decay rate

# Basis
b_cav = FockBasis(4)
b_mech = FockBasis(10)

# Operators Cavity
a = destroy(b_cav) ⊗ one(b_mech)
at = create(b_cav) ⊗ one(b_mech)
an = number(b_cav) ⊗ one(b_mech)

# Operators Oscillator
b = one(b_cav) ⊗ destroy(b_mech)
bt = one(b_cav) ⊗ create(b_mech)
bn = one(b_cav) ⊗ number(b_mech)

function H(
	ω_mech=ω_mech, Δ=Δ, g=g, η=η, κ=κ, γ=γ,
	b_cav=b_cav, b_mech=b_mech,
	a=a, at=at, b=b, bt=bt
)
	H_cav = -Δ * at * a + η * (a + at)
	H_mech = ω_mech * bt * b + γ * (bt * b)
	H_int = g * (at * a) * (bt + b)
	return H_cav + H_mech + H_int	
end

ψ0 = fockstate(b_cav,1) ⊗ fockstate(b_mech,1)

# Time evolution
T = [0:0.2:20;]
@time tout, ψt = timeevolution.master(T, ψ0, H(), [a]; rates=[1])

f, ax, p = lines(T, real(expect(an, ψt)), label="optical")
lines!(T, real(expect(bn, ψt)), label="mechanical")
ax.xlabel="Time"
ax.ylabel="<n>"
ylims!(low=0)
xlims!(low=0)
axislegend(position=:rc)
f