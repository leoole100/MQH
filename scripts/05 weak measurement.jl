include("00 functions.jl")
using QuantumOptics, CairoMakie, LinearAlgebra, StatsBase

# Define the Hilbert space (two-level system)
basis = SpinBasis(1//2)
κ=0.5


function measure_weak(ρ::Operator, M::Operator, p::Float64)
	"""Performs a weak measurement on the state ρ using observable M with strength p.
		 When p=1, the measurement collapses the state to an eigenstate of M.
	"""
	
	# Compute eigenvalues and eigenvectors of M
	eigenvals, eigenvecs = eigen(dense(M).data)
	
	# Get Hilbert space basis
	b = ρ.basis_l
	
	# Define the "no measurement" Kraus operator
	M0 = sqrt(1 - p) * identityoperator(b)

	# Define the weak measurement Kraus operators
	Kraus_ops = [sqrt(p) * projector(Ket(b, eigenvecs[:, i])) for i in 1:length(eigenvals)]
	
	# Compute probabilities
	probs = [real(tr(K * ρ * dagger(K))) for K in Kraus_ops]
	p0 = real(tr(M0 * ρ * dagger(M0)))  # Probability of no measurement
	push!(probs, p0)  # Include no-measurement probability

	println(probs)

	# Choose measurement outcome randomly
	outcome_index = sample(1:length(probs), Weights(probs))
	
	# Apply the corresponding Kraus operator
	if outcome_index == length(probs)  # No significant update (M0 case)
		ρ_new = M0 * ρ * dagger(M0) / tr(M0 * ρ * dagger(M0))
		m = sum(eigenvals .* probs[1:end-1])  # Expectation value
	else  # Collapse toward an eigenstate
			K_chosen = Kraus_ops[outcome_index]
			ρ_proj = K_chosen * ρ * dagger(K_chosen) / tr(K_chosen * ρ * dagger(K_chosen))
			ρ_new = (1 - p) * ρ + p * ρ_proj
			ρ_new = (ρ_new + dagger(ρ_new)) / 2

			# Generate noisy measurement outcome
			noise = randn() * sqrt(1-p)
			m = eigenvals[outcome_index] + noise
	end

	return m, ρ_new
end

# Define a Hilbert space
b = SpinBasis(1//2)

# Define an observable with multiple eigenvalues
σz = sigmaz(b)  # Pauli Z matrix (eigenvalues: +1, -1)

# Initial state: Superposition (|0⟩ + |1⟩) / sqrt(2)
ψ0 = normalize(spinup(b) + spindown(b))
ρ = dm(ψ0)

# Weak measurement
p = 0.5  # Weak measurement strength
m, ρ_new = measure_weak(ρ, σz, p)

println("Measurement outcome: ", m)
println("Updated state: ", ρ_new)

# %%
function evolution(
	times, meas_times, p;
	H = dense(sigmaz(basis)),
	ψ0 = spinup(basis),
	C = sigmaz(basis),
	J = κ*sigmam(basis)
)
	time = 0
	state = dm(ψ0)

	tout = Float64[]
	states = []

	sort!(meas_times)
	dt = times[2]

	for i in 1:length(meas_times)
		evo = timeevolution.master(time:dt:meas_times[i], state, H, [J])
		append!(tout, evo[1])
		append!(states, evo[2])
		time = tout[end]+dt
		# _, state = measure(C, Ket(basis, diag(states[end].data)))
		_, state = measure_weak(C, states[end], p[i])
	end


	evo = timeevolution.master(time:dt:maximum(times), state, H, [J])
	append!(tout, evo[1])
	append!(states, evo[2])
	
	return tout, states
end


# evo = evolution(0:0.01:10, [], [])
# strong = evolution(0:0.01:10, [1, 4], [1., 1.])
weak = evolution(0:0.01:10, [1, 4], [0.1,0.1])


f = Figure(size=(fullsize[1]*0.8, fullsize[2]))
a = Axis(f[1,1], 
	xlabel="time",
	yticks=([-1,1, 0], ["↓", "↑", ""]),
)
plot_sim(evo, Makie.wong_colors()[1])
plot_sim(strong, Makie.wong_colors()[2])
plot_sim(weak, Makie.wong_colors()[3])
hidexdecorations!(a)
# save("../figures/05 discreate.pdf", f)
display(f)

# %%