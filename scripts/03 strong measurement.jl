include("00 functions.jl")
using QuantumOptics, CairoMakie, LinearAlgebra, StatsBase

# Define the Hilbert space (two-level system)
basis = SpinBasis(1//2)

function evolve_and_measure(time_array, measurement_time, H, C, ψ0)
	time_pre = 0:time_array[2]:(measurement_time - time_array[1])
	pre = timeevolution.schroedinger(time_pre, ψ0, H)
	
	ψ = pre[2][end]
	outcome, ψ = measure(C, ψ)
	
	time_post = measurement_time:time_array[2]:maximum(time_array)
	post = timeevolution.schroedinger(time_post, ψ, H)
	
	# Combine the results (evolution before and after measurement)
	full_time = vcat(pre[1], post[1])
	full_state = vcat(pre[2], post[2])
	
	return full_time, full_state
end

function free_and_measured(time_array, measure_time; H, C, ψ0)
	free = timeevolution.schroedinger(time_array, ψ0, H)
	measured = evolve_and_measure(time_array, measure_time, H, C, ψ0)
	return free, measured
end

sims = free_and_measured(0:0.01:2, 1,
	H = dense(sigmaz(basis)),
	C = sigmaz(basis),
	ψ0 = normalize(spinup(basis) + spindown(basis))
)

f = Figure(size=halfsize)
a = Axis(f[1,1], 
	xlabel="time",
	yticks=([-1,1, 0], ["↓", "↑", ""]),
)
plot_sim(sims[1], Makie.wong_colors()[1])
hidexdecorations!(a)
save("../figures/03 free.pdf", f)
display(f)

f = Figure(size=halfsize)
a = Axis(f[1,1], 
	xlabel="time",
	yticks=([-1,1, 0], ["↓", "↑", ""]),
)
plot_sim(sims[2], Makie.wong_colors()[2])
hidexdecorations!(a)
save("../figures/03 measurement.pdf", f)
display(f)