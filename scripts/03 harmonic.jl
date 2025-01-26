include("00 functions.jl")
using QuantumOptics, CairoMakie

# %%
b = FockBasis(4)
a = destroy(b)
ad = create(b)

H_harmonic = ad*a

C_harmonic = ad+a

ψ0 = coherentstate(b, 1)

timespan = 2*π
free = timeevolution.schroedinger([0:0.1:timespan;], ψ0, H_harmonic)
measured = timeevolution.schroedinger([0:0.1:timespan;], normalize(C_harmonic*ψ0), H_harmonic)


function plot_sim(sim, O=proj)
	expectation = expect(O, sim[2])
	uncertainty = sqrt.(expect(O*O, sim[2]) .- expectation.^2)
	band!(sim[1],
		real(expectation-uncertainty),
		real(expectation+uncertainty),
		alpha=.1
	)
	return lines!(sim[1], real(expectation))
end

f = Figure()
ax = Axis(f[1,1])
plot_sim(free, ad+a)
plot_sim(measured, ad+a)
hidexdecorations!(ax, grid=false)

ax = Axis(f[2, 1])
plot_sim(free, 1im*(ad-a))
plot_sim(measured, 1im*(ad-a))

rowgap!(f.layout, 0)
f