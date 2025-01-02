using QuantumOptics, CairoMakie, LinearAlgebra, FFTW, Random
include("00_functions.jl")
cd(@__DIR__)

x = δa + dagger(δa)
ψ0 = tensor(fockstate(optical_space, 0), fockstate(mechanical_space, 0))

function calc_S_xx(HJC)
	times_full = 0:0.25:20

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

function S_xx_res(HJC)
	frequencies, S_xx = calc_S_xx(HJC)
	return S_xx[findmin(x->abs(x-1), frequencies)[2]]	
end

f = Figure(size=fullsize)
a = Axis(f[1,1], xlabel="frequency in ωₘ", ylabel="Sₓₓ")

drive_amplitude = range(0, .5, 5)
@time for da in drive_amplitude
	lines!(calc_S_xx(HJC(drive_amplitude=da.*ω_m))...,#%%
	function power_sweep(n_th)
		drive_amplitude = logrange(1e-10, .5, 10)
		S = [S_xx_res(HJ(drive_amplitude=a.*ω_m, n_th=n_th, γ_m=0.1*ω_m)) for a in drive_amplitude]
		return drive_amplitude, S
	end
	
	n_th = [2]
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
	save("../figures/03 power sweep.pdf", f)
	f
		color=da, colorrange=extrema(drive_amplitude),
		label="$(round(da, digits=2))" 
	)
end
axislegend("Drive Strength")
ylims!(a, low=0)
xlims!(a, low=0)
save("../figures/04 power spectrum.pdf", f)
f

#%%
function power_sweep(n_th)
	drive_amplitude = logrange(1e-10, .5, 10)
	S = [S_xx_res(HJC(drive_amplitude=a.*ω_m, n_th=n_th, γ_m=0.1*ω_m)) for a in drive_amplitude]
	return drive_amplitude, S
end

n_th = [2]
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