using QuantumOptics, CairoMakie, LinearAlgebra, FFTW, Random
include("00_functions.jl")
cd(@__DIR__)

x = δa + dagger(δa)
ψ0 = tensor(fockstate(optical_space, 0), fockstate(mechanical_space, 0))

function calc_S_xx(HJC; time=20)
	times_full = 0:0.01:time
	
	tout, result_full = stochastic.master(times_full, ψ0, HJC..., dt=1e-3)
	
	cutoff = 0
	result = result_full[times_full .>= cutoff]
	times = times_full[times_full .>= cutoff]
	S_xx = abs2.(fft(expect(x, result_full))) 
	
	frequencies = fftfreq(length(times), 1/(times[2]-times[1]))
	S_xx = S_xx[frequencies .> 0]
	frequencies = frequencies[frequencies .> 0]
	
	return frequencies, S_xx
end

function S_xx_res(HJC; time=20)
	frequencies, S_xx = calc_S_xx(HJC, time=time)
	return S_xx[findmin(x->abs(x-1), frequencies)[2]]	
end

# %%
times = 0:0.01:20

c = HJC(
	ω_m = 2π * 1,      					# Mechanical frequency (Hz)
	κ = 0.1 * ω_m,        			# Cavity linewidth (Hz)
	g0 = 0.1*ω_m,        				# Single-photon coupling (Hz)
	Δ = ω_m,            	 			# Detuning (drive on red sideband)
	drive_amplitude = 0.2*ω_m,  # Drive amplitude (Hz)
	γ_m = 0.05 * ω_m,      		# Mechanical damping rate (Hz)
	n_th = 0,           				# Thermal occupation number
	η=.1
)

c = HJC(
	ω_m = 2π * 1,      					# Mechanical frequency (Hz)
	κ = 1.0 * ω_m,        			# Cavity linewidth (Hz)
	g0 = 0.05*ω_m,        				# Single-photon coupling (Hz)
	Δ = ω_m,            	 			# Detuning (drive on red sideband)
	drive_amplitude = 0.01*ω_m,  # Drive amplitude (Hz)
	γ_m = 0.01 * ω_m,      		# Mechanical damping rate (Hz)
	n_th = 0.1,           				# Thermal occupation number
	η=1
)
	
tout, results = timeevolution.master(times, ψ0, c[1:2]...)
tout, results_stoch = stochastic.master(times, ψ0, c..., dt=.001)

f = Figure(size=fullsize)
a = Axis(f[1,1], ylabel="<x>", xlabel="T")
lines!(times, abs2.(expect(x, results)), label="free system")
lines!(times, abs2.(expect(x, results_stoch)), label="measurement")
axislegend()
save("../figures/04 time evolution.pdf", f)
f


#%%
f = Figure(size=fullsize)
a = Axis(f[1,1], xlabel="frequency in ωₘ", ylabel="Sₓₓ")

drive_amplitude = range(0, .5, 3)
@time for da in drive_amplitude
	lines!(calc_S_xx(HJC(drive_amplitude=da.*ω_m))...,
		color=da, colorrange=extrema(drive_amplitude),
		label="$(round(da, digits=2))" 
	)
end
axislegend("Drive Strength")
ylims!(a, low=0)
xlims!(a, 0, 2)
save("../figures/04 power spectrum.pdf", f)
f

#%%
function power_sweep(n_th)
	drive_amplitude = logrange(1e-5, 10, 10)
	S = [S_xx_res(HJC(drive_amplitude=a.*ω_m, n_th=n_th, η=1); time=20) for a in drive_amplitude]
	return drive_amplitude, S
end

n_th = [0]
@time res = power_sweep.(n_th);

f = Figure(size=fullsize)
a = Axis(f[1,1], 
	xlabel="laser Power in ωₘ",
	ylabel="Sₓₓ(ωₘ)",
	xscale=log10, yscale=log10
)
# for (n, r) in zip(n_th, res)
# 	scatterlines!(r..., color=n, colorrange=extrema(n_th), label="n=$(n)")
# end
# axislegend(position=:lt)
scatterlines!(res[1][1], abs.(res[1][2]))
save("../figures/04 power sweep.pdf", f)
f