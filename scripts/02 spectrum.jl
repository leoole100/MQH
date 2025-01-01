using QuantumOptics, CairoMakie, LinearAlgebra, FFTW
include("00_functions.jl")
cd(@__DIR__)

#%%
optical_space = FockBasis(10)
mechanical_space = FockBasis(4)
hilbert_space = optical_space ⊗ mechanical_space

δa = destroy(optical_space) ⊗ one(mechanical_space)
δb = one(optical_space) ⊗ destroy(mechanical_space)

x = δa + dagger(δa)

ψ0 = tensor(fockstate(optical_space, 0), fockstate(mechanical_space, 0))

function HJ(;
	ω_m = 2π * 1,      					# Mechanical frequency (Hz)
	κ = 0.1 * ω_m,        			# Cavity linewidth (Hz)
	g0 = 0.1*ω_m,        				# Single-photon coupling (Hz)
	Δ = ω_m,            	 			# Detuning (drive on red sideband)
	drive_amplitude = 0.2*ω_m,  # Drive amplitude (Hz)
	γ_m = 0.05 * ω_m,      		# Mechanical damping rate (Hz)
	n_th = 2,           				# Thermal occupation number
)

	# Steady-state amplitudes (classical solution)
	α = drive_amplitude / (κ/2 - 1im * Δ)  # Steady-state optical field

	# Effective coupling
	G = g0 * abs(α)  # Linearized coupling strength
	
	# Linearized Hamiltonian
	H_c = Δ * dagger(δa) * δa
	H_m = ω_m * dagger(δb) * δb
	H_int = -G * (δa + dagger(δa)) * (δb + dagger(δb))
	H_pump = drive_amplitude * (δa + dagger(δa))
	H = H_c + H_m + H_int + H_pump

	# Collapse operators
	C_optical = √κ * δa
	C_mechanical = √(γ_m * (n_th + 1)) * δb
	C_thermal = √(γ_m * n_th) * dagger(δb)
	J = [C_optical, C_mechanical, C_thermal]

	return H, J
end

times_full = 0:0.25:50
@time tout, result_full = timeevolution.master(times_full, ψ0, HJ()...);

cutoff = 0
result = result_full[times_full .>= cutoff]
times = times_full[times_full .>= cutoff]
S_xx = abs2.(fft(expect(x, result)))
frequencies = fftfreq(length(times), 1/(times[2]-times[1]))
S_xx = S_xx[frequencies .> 0]
frequencies = frequencies[frequencies .> 0]

f = Figure(size=fullsize)
at = Axis(f[1,1], xlabel="Time",  ylabel="x")
af = Axis(f[2,1], xlabel="Frequency in ω_m", ylabel="Sₓₓ")
ylims!(af, low=0)
xlims!(af, low=0 )

lines!(at, times, abs2.(expect(x, result)))
lines!(af, frequencies, S_xx)
f