include("00 functions.jl")
using QuantumOptics, CairoMakie, LinearAlgebra, StatsBase, ColorSchemes

# Define the Hilbert space (two-level system)
basis = SpinBasis(1//2)

function evolution(;
	times=0:0.01:10,
	H = dense(sigmaz(basis)),
	ψ0 = spinup(basis),
	eta=0.4,
	κ=0.5,
	C = eta*sigmaz(basis),
	J = κ*sigmam(basis),
	dt=1e-3
)

	evo = stochastic.master(times, dm(ψ0), H, [J], [C]; dt=dt)
	meas = expect(C./eta, evo[2]) .+ randn(Complex{Float64}, length(evo[2])) .* sqrt(dt/eta)

	return evo, meas
end

function get_value(measured; O=sigmaz(basis))
	return real.(expect_uncert(O, measured[2]))[1]
end

eta = [0.01, 0.1, 1]
evolutions = [evolution(eta=e) for e in eta]

# %%

f = Figure(size=(0.8*fullsize[1], 0.9*fullsize[2]))
a = Axis(f[1,1], 
	xlabel="time",
	yticks=([-1,1, 0], ["↓", "↑", ""]),
	ylabel="Evolution"
)
b = Axis(f[2,1],
	yticks=([-1,1, 0], ["↓", "↑", ""]),
	ylabel="Measured"
) 

for (e, evolution) in zip(eta, evolutions)
	evo = evolution[1]
	lines!(a, evo[1], get_value(evo), label="$e", color=log10(e), colorrange=extrema(log10.(eta)))
	lines!(b, evo[1], real.(evolution[2]), label="$e", color=log10(e), colorrange=extrema(log10.(eta)))
end


Colorbar(f[:,2], colorrange=extrema(log10.(eta)), label="log η", ticks=log10.(eta))
hidexdecorations!(a)
hidexdecorations!(b)
rowgap!(f.layout, 0)
colgap!(f.layout, 0)
save("../figures/06 sme.pdf", f)
display(f)

