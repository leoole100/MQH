# make figure style
using Makie
inch = 96
pt = 4/3
cm = inch/2.54
update_theme!(
	fontsize=10pt,
	px_per_unit = 10,
	fonts = (;
		regular="Roboto Regular",
		bold="Roboto Bold",
	),
)
halfsize = (2inch, (2/1.33)inch)
fullsize = (4inch, (4/1.33)inch)


# Define the basis
using QuantumOptics
optical_space = FockBasis(10)
mechanical_space = FockBasis(4)
hilbert_space = optical_space ⊗ mechanical_space

# and the commonly used operators
δa = destroy(optical_space) ⊗ one(mechanical_space)
δb = one(optical_space) ⊗ destroy(mechanical_space)

# define H and J in terms of Parameters
ω_m = 2π * 1
function HJ(;
	ω_m = 2π * 1,      					# Mechanical frequency (Hz)
	κ = 0.1 * ω_m,        			# Cavity linewidth (Hz)
	g0 = 0.1*ω_m,        				# Single-photon coupling (Hz)
	Δ = ω_m,            	 			# Detuning (drive on red sideband)
	drive_amplitude = 0.2*ω_m,  # Drive amplitude (Hz)
	γ_m = 0.05 * ω_m,      		# Mechanical damping rate (Hz)
	n_th = 2,           				# Thermal occupation number
	δa = δa,
	δb = δb,
)

	# Steady-state amplitudes (classical solution)
	α = drive_amplitude / (κ/2 - 1im * Δ)  # Steady-state optical field

	# Effective coupling
	G = g0 * abs(α)  # Linearized coupling strength
	
	# Linearized Hamiltonian
	H_c = Δ * dagger(δa) * δa
	H_m = ω_m * dagger(δb) * δb
	H_int = -G * (δa + dagger(δa)) * (δb + dagger(δb))
	H_pump = im*drive_amplitude * (dagger(δa) - δa)
	H = H_c + H_m + H_int + H_pump

	# Collapse operators
	C_optical = √κ * δa
	C_mechanical = √(γ_m * (n_th + 1)) * δb
	C_thermal = √(γ_m * n_th) * dagger(δb)
	J = [C_optical, C_mechanical, C_thermal]

	return H, J
end

function HJC(;
	ω_m = 2π * 1,      					# Mechanical frequency (Hz)
	κ = 0.1 * ω_m,        			# Cavity linewidth (Hz)
	g0 = 0.1*ω_m,        				# Single-photon coupling (Hz)
	Δ = ω_m,            	 			# Detuning (drive on red sideband)
	drive_amplitude = 0.2*ω_m,  # Drive amplitude (Hz)
	γ_m = 0.05 * ω_m,      		# Mechanical damping rate (Hz)
	n_th = 2,           				# Thermal occupation number
)

	H, J = HJ(
		ω_m = ω_m,
		κ = κ,
		g0 = g0,
		Δ = Δ,
		drive_amplitude = drive_amplitude,
		γ_m = γ_m,
		n_th = n_th
	)

	C = sqrt(κ)*(δa+dagger(δa))

	return H, J, [C]
end