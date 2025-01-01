using QuantumOptics, CairoMakie, LinearAlgebra
include("00_functions.jl")
cd(@__DIR__)

# Parameters
ω_m = 2π * 1      # Mechanical frequency (Hz)
κ = 0.1 * ω_m        # Cavity linewidth (Hz)
γ_m = 0.001 * ω_m      # Mechanical damping rate (Hz)
g0 = 0.01*ω_m        # Single-photon coupling (Hz)
Δ = ω_m             # Detuning (drive on red sideband)
drive_amplitude = 100*ω_m  # Drive amplitude (Hz)
n_th = 0           # Thermal occupation number


# Steady-state amplitudes (classical solution)
α = drive_amplitude / (κ/2 - 1im * Δ)  # Steady-state optical field

# Effective coupling
G = g0 * abs(α)  # Linearized coupling strength

# Operators
optical_space = FockBasis(10)
mechanical_space = FockBasis(10)
hilbert_space = optical_space ⊗ mechanical_space

δa = destroy(optical_space) ⊗ one(mechanical_space)
δb = one(optical_space) ⊗ destroy(mechanical_space)

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
collapse_ops = [C_optical, C_mechanical, C_thermal]

# Initial state: vacuum for fluctuations
ψ0 = tensor(fockstate(optical_space, 0), fockstate(mechanical_space, 0))

# Time evolution
times = 0:0.01:3
@time tout, result = timeevolution.master(times, ψ0, H, collapse_ops);

# %%
# look at timeevolution
function numberStates(state, index::Integer)
	trace = ptrace.(state, index)
	return hcat([diag(normalize(t).data) for t in trace]...)'	
end

N_opt = [expect(p, dagger(δa)*δa) for p in result]
N_mech = [expect(p, dagger(δb)*δb) for p in result]

f = Figure(size=fullsize)
a = Axis(f[1,1], xlabel="Time in 2π/ωₘ", ylabel="n")

lines!(times, abs.(N_opt), label="Optical")
heatmap!(times, 0:optical_space.N, abs.(numberStates(result, 2)), 
 colormap=[:transparent, Makie.wong_colors()[1]]
)

lines!(times, abs.(N_mech), label="Mechanical")

axislegend(position=:rb, padding=0, framevisible=false)

save("../figures/02 time evolution.pdf", f)
f
