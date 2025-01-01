using QuantumOptics, CairoMakie, LinearAlgebra
include("00_functions.jl")
cd(@__DIR__)

# Parameters
ω_m = 2π * 1
H, J = HJ(
	ω_m = ω_m,      					# Mechanical frequency (Hz)
	κ = 0.1 * ω_m,        			# Cavity linewidth (Hz)
	g0 = 0.1*ω_m,        				# Single-photon coupling (Hz)
	Δ = ω_m,            	 			# Detuning (drive on red sideband)
	drive_amplitude = 0.2*ω_m,  # Drive amplitude (Hz)
	γ_m = 0.003 * ω_m,      			# Mechanical damping rate (Hz)
	n_th = 2,           				# Thermal occupation number
)

# Initial state: vacuum for fluctuations
ψ0 = tensor(fockstate(optical_space, 0), fockstate(mechanical_space, 1))

# Time evolution
times = 0:0.2:50
@time tout, result = timeevolution.master(times, ψ0, H, J);

# look at timeevolution
function numberStates(state, index::Integer)
	trace = ptrace.(state, index)
	return hcat([diag(normalize(t).data) for t in trace]...)'	
end

f = Figure(size=fullsize)
aO = Axis(f[1,1], ylabel="Optical")
aM = Axis(f[2,1], xlabel="Time in 2π/ωₘ", ylabel="Mechanical")
hidexdecorations!(aO, grid=false)
linkxaxes!(aM, aO)
ylims!(aO, low=0)
ylims!(aM, low=0)

lines!(aO, times, abs.(expect(dagger(δa)*δa, result)))
lines!(aM, times, abs.(expect(dagger(δb)*δb, result)))
#heatmap!(aO, times, 0:optical_space.N, abs.(numberStates(result, 2)), 
#  colormap=[:transparent, Makie.wong_colors()[1]]
# )

save("../figures/01 cooling.pdf", f)
f