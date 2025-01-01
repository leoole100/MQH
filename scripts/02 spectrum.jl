using QuantumOptics, CairoMakie, LinearAlgebra, FFTW
include("00_functions.jl")
cd(@__DIR__)

x = δa + dagger(δa)
ψ0 = tensor(fockstate(optical_space, 0), fockstate(mechanical_space, 0))

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
af = Axis(f[2,1], xlabel="Frequency in ωₘ", ylabel="Sₓₓ")
ylims!(af, low=0)
xlims!(af, low=0 )

lines!(at, times, abs2.(expect(x, result)))
lines!(af, frequencies, S_xx)
save("../figures/02 spectrum.pdf", f)
f