include("00 functions.jl")
using QuantumOptics, CairoMakie, LinearAlgebra, StatsBase

# Define the Hilbert space (two-level system)
basis = SpinBasis(1//2) ⊗ FockBasis(3)

a = destroy(bc) ⊗ one(ba)
ad = create(bc) ⊗ one(ba)
σ⁺ = one(bc) ⊗ sigmap(ba)
σ = one(bc) ⊗ sigmam(ba)


function evolution(
	times, meas_times;
	κ=0.5,
	H = dense(sigmaz(basis)),
	ψ0 = spinup(basis) ⊗ fockstate(bc, 0),
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


evo = evolution(0:0.01:10, [])
meas = evolution(0:0.01:10, [1, 6])


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

times = 0:0.1:10
evo = evolution(times, [])
meas = evolution(times, collect(1:1:9))

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