using QuantumOptics, CairoMakie, LinearAlgebra, FFTW
include("00_functions.jl")
cd(@__DIR__)

x = δa + dagger(δa)
ψ0 = tensor(fockstate(optical_space, 0), fockstate(mechanical_space, 0))


function S_xx_res(HJ)
	times_full = 0:0.25:20
	tout, result_full = timeevolution.master(times_full, ψ0, HJ...);

	cutoff = 0
	result = result_full[times_full .>= cutoff]
	times = times_full[times_full .>= cutoff]
	S_xx = abs2.(fft(expect(x, result)))
	frequencies = fftfreq(length(times), 1/(times[2]-times[1]))
	S_xx = S_xx[frequencies .> 0]
	frequencies = frequencies[frequencies .> 0]

	return S_xx[findmin(x->abs(x-1), frequencies)[2]]
end


drive_amplitude = logrange(0.0001, .5, 10)
@time S = [S_xx_res(HJ(drive_amplitude=a.*ω_m)) for a in drive_amplitude]

f = Figure(size=fullsize)
a = Axis(f[1,1], 
	xlabel="laser Power in ωₘ",
	ylabel="Sₓₓ(ωₘ)",
	xscale=log10, yscale=log10
)
scatterlines!(drive_amplitude, S)
save("../figures/03 power sweep.pdf", f)
f