include("00 functions.jl")
using QuantumOptics, CairoMakie, LinearAlgebra, StatsBase

# Define the Hilbert space (two-level system)
basis = SpinBasis(1//2)

# Define the system's Hamiltonian (e.g., a simple Rabi oscillation)
H = dense(sigmax(basis))

# Initial state (e.g., |ψ⟩ = |↑⟩)
ψ0 = spinup(basis)

# measurement operator
C = sigmaz(basis)


time = 0
ψ = ψ0
times = []
states = []
total_time = 1π
N = 5*total_time/π
for i in 1:N
	t = time:0.01:time+total_time/N
	time = t[end]+0.01
	tout, s = timeevolution.schroedinger(t, ψ, H)
	outcome, ψ = measure(C, ψ)

	append!(times, tout)
	append!(states, s)
end

f = Figure(size=fullsize)
a = Axis(f[1,1], 
	xlabel="time", ylabel=L"\sigma_z",
	yticks=([-1,1, 0], ["↓", "↑", ""]),
)
O = C 

obs = real.(expect_uncert(O, states))
lines!(times, obs[1])
band!(float.(times), obs[1]-obs[2], obs[1]+obs[2], alpha=.2)

vlines!((total_time/N)*(1:N))

xlims!(a, extrema(times))
ylims!(a, low=-1)
save("../figures/04 discreate zeno.pdf", f)
display(f)