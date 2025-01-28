#==
time evolution of the rabi oscillation with a single, discrete measurement
==#
include("00 functions.jl")
using CairoMakie, QuantumOptics

# Parameters
Nc = 10
κ = 1.0
g = .2κ
β = .1*sqrt(κ)
η = .9
gs = 0.1κ

# Basis and operators
bc = FockBasis(Nc)
ba = SpinBasis(1//2)

a = destroy(bc) ⊗ one(ba)
ad = create(bc) ⊗ one(ba)
σ⁺ = one(bc) ⊗ sigmap(ba)
σ = one(bc) ⊗ sigmam(ba)

H(;g=g, gs=gs, β=β) = g*ad*a*σ⁺*σ + gs*(σ⁺+σ) - -1.0im*β*(ad-a)
J(;κ=κ) = [κ*a]
C(;κ=κ, η=η) = [sqrt(κ*η)*a]

ψ0 = coherentstate(bc, β) ⊗ spindown(ba)

# make time evolution
time_end = 2*π/gs
free = timeevolution.master([0:0.1:time_end;], ψ0, H(), J())

time_measurement = 2
measured_start = timeevolution.master([0:0.1:time_measurement;], ψ0, H(), J())
measured = C()[1] * measured_start[2][end]
measured = normalize(measured)
measured_after = timeevolution.master([time_measurement:0.1:time_end;], measured, H(), J())

@time continous = stochastic.master([0:0.1:time_end;], ψ0, H(), J(), C(), dt=1e-3);

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

#%%
f = Figure(size=fullsize)

ax = Axis(f[1,1], ylabel="<σ⁺σ>")
proj =  σ⁺*σ
plot_sim(free)
plot_sim(measured_after)
plot_sim(continous)
hidexdecorations!(ax, grid=false)

ax = Axis(f[2,1], xlabel="time", ylabel="<σ⁺ + σ>")
proj =  σ⁺ + σ
plot_sim(free)
plot_sim(measured_after)
plot_sim(continous)
hidexdecorations!(ax, grid=false)

ax = Axis(f[3,1], xlabel="time", ylabel="<σ⁺ - σ>")
proj =  σ⁺ - σ
axislegend(ax, [
	plot_sim(free),
	plot_sim(measured_after),
	plot_sim(continous)
],
["free", "single", "continous"],
margin=[0,0,0,0], padding=2)

rowgap!(f.layout, 0)
save("../figures/02 rabi w measurement.pdf", f)
f