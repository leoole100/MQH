#==
Just look at the cavity without any noise effects.
==#

using QuantumOptics
using CairoMakie

#%%
G = 1
Δ = 0
Ω = 1 

basis_a = FockBasis(20)
basis_b = FockBasis(50)

a = destroy(basis_a) ⊗ one(basis_b)
b = one(basis_a) ⊗ destroy(basis_b)

x = b + dagger(b)
na = dagger(a)*a
nb = dagger(b)*b

H = -Δ*dagger(a)*a - G*x*dagger(a)*a + Ω*dagger(b)*b

ψ0 = coherentstate(basis_a, 3) ⊗ fockstate(basis_b, 0)

# display the state probabilities <i,j|ψ>
function selector(basis::Basis, n::Integer)
	z = zeros(length(basis))
	z[n] = 1 
	return Bra(basis, z)
end

function state_prob(ψ::Ket, basis_a=basis_a, basis_b=basis_b)
	return [
		abs(selector(basis_a, i) ⊗ selector(basis_b, j) * ψ)^2
		for i in 1:length(basis_a), j in 1:length(basis_b)]
end

function plot_state(state_prob::Matrix)
	f,a,p = heatmap(0:length(basis_a)-1, 0:length(basis_b)-1, state_prob)
	a.xlabel = "a"
	a.ylabel = "b"
	f
end
plot_state(s::Ket, basis_a=basis_a, basis_b=basis_b) = plot_state(state_prob(s, basis_a, basis_b))

# plot_state(ψ0)


# Time evolution
T = [0:0.1:10;]
@time tout, ψt = timeevolution.schroedinger(T, ψ0, H)

f = Figure()
at = Axis(f[1, 1:2], xlabel="time", ylabel="<n>")

# lines!(T, abs.(expect(na, ψt)), label="optical")
# lines!(T, abs.(expect(nb, ψt)), label="mechanical")
lines!(T, abs.(expect(H, ψt)), label="H")
ylims!(low=0)
# xlims!(low=0)
axislegend(position=:rc)

a0 = Axis(f[2,1], title="start", ylabel="b", xlabel="a")
heatmap!(a0, 0:length(basis_a)-1, 0:length(basis_b)-1, state_prob(ψt[1]))
ae = Axis(f[2,2], title="end", xlabel="a")
heatmap!(ae, 0:length(basis_a)-1, 0:length(basis_b)-1, state_prob(ψt[end]))
f