include("00 functions.jl")
using QuantumOptics, CairoMakie, LinearAlgebra, StatsBase

# %%
# Define the Hilbert space (two-level system)
basis = SpinBasis(1//2)

# Define the system's Hamiltonian (e.g., a simple Rabi oscillation)
H = dense(sigmax(basis))

# Initial state (e.g., |ψ⟩ = |↑⟩)
ψ0 = spinup(basis)

# measurement operator
C = sigmaz(basis)

function measure(C::Operator, ψ::Ket)
	eigenvals, eigenvecs = eigen(dense(C).data)
	
	probs = [abs2(dagger(ψ)*Ket(basis,eigenvecs[:, i])) for i in 1:length(eigenvals)]
	
	outcome_index = sample( 1:length(eigenvals), Weights(probs))
	outcome_value = eigenvals[outcome_index]
	ψ_post = normalize(Ket(basis, eigenvecs[:, outcome_index]))
	
	return outcome_value, ψ_post
end

function evolve_and_measure(time_array, measurement_time)
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

times = 0:0.01:π
measure_time =  π/4
measured = evolve_and_measure(times, measure_time)
free = timeevolution.schroedinger(times, ψ0, H)

function expect_uncert(O::Operator, ψs)
	expectation = expect(O, ψs)
	uncertainty = sqrt.(expect(O*O, ψs) .- expectation.^2)
	return expectation, uncertainty
end

f = Figure(size=fullsize)
a = Axis(f[1,1], 
	xlabel="time", ylabel=L"\sigma_z",
	yticks=([-1,1], ["↓", "↑"])
)
O = C 

obs = real.(expect_uncert(O, measured[2]))
lines!(measured[1], obs[1], label="measured")
band!(measured[1], obs[1]-obs[2], obs[1]+obs[2], alpha=.2)

obs = real.(expect_uncert(O, free[2]))
lines!(free[1], obs[1], label="free")
band!(free[1], obs[1]-obs[2], obs[1]+obs[2], alpha=.2)

# vlines!(measure_time)
axislegend(framevisible=false, position=:rb)
xlims!(a, extrema(times))
save("../figures/03 measurement.pdf", f)
display(f)