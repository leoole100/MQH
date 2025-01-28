include("00 functions.jl")
using QuantumOptics, CairoMakie, LinearAlgebra, StatsBase

# Define the Hilbert space (two-level system)
basis = SpinBasis(1//2)
κ=0.5

H = dense(sigmaz(basis))
ψ0 = spinup(basis)
C = sigmaz(basis)
J = κ*sigmam(basis)

function evolution(
	times, meas_times;
	H = dense(sigmaz(basis)),
	ψ0 = spinup(basis),
	C = sigmaz(basis),
	J = κ*sigmam(basis)
)
	time = 0
	state = ψ0

	tout = Float64[]
	states = []

	sort!(meas_times)
	dt = times[2]

	for i in 1:length(meas_times)
		evo = timeevolution.master(time:dt:meas_times[i], state, H, [J])
		append!(tout, evo[1])
		append!(states, evo[2])
		time = tout[end]+dt
		_, state = measure(C, Ket(basis, diag(states[end].data)))
	end


	evo = timeevolution.master(time:dt:maximum(times), state, H, [J])
	append!(tout, evo[1])
	append!(states, evo[2])
	
	return tout, states
end


evo = evolution(0:0.1:10, [])
meas = evolution(0:0.1:10, [2, 5])


f = Figure(size=(fullsize[1]*0.8, fullsize[2]))
a = Axis(f[1,1], 
	xlabel="time",
	yticks=([-1,1, 0], ["↓", "↑", ""]),
)
plot_sim(evo, Makie.wong_colors()[1])
plot_sim(meas, Makie.wong_colors()[2])
hidexdecorations!(a)
save("../figures/04 discreate.pdf", f)
display(f)


# %%

evo = evolution(0:0.1:10, [])
meas = evolution(0:0.1:10, collect(1:2:9))

f = Figure(size=(fullsize[1], fullsize[2]))
a = Axis(f[1,1], 
	xlabel="time",
	yticks=([-1,1, 0], ["↓", "↑", ""]),
)
plot_sim(evo, Makie.wong_colors()[1])
plot_sim(meas, Makie.wong_colors()[2])
hidexdecorations!(a)
save("../figures/04 discreate zeno.pdf", f)
display(f)