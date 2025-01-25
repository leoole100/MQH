include("00 functions.jl")

# %%

ψ0 = tensor(coherentstate(optical_space, 0), coherentstate(mechanical_space, 0))

# plot the initial state
function plot_wigner(ρ, index, s=Inf)
	t = ptrace(ρ, index)
	s = minimum([s, (size(t, 1)-1)])
	x = -s:s
	y = -s:s
	return x, y, wigner(t, x, y)
end

function plot_wigner(ρ)
	f = Figure()
	a = Axis(f[1,1], title="Mechanical", aspect=DataAspect(), xlabel="Q", ylabel="P")
	heatmap!(plot_wigner(ρ, 1, 5)..., colormap=:diverging, colorrange=(-.01,.01))
	a = Axis(f[1,2], title="Optical", aspect=DataAspect(), xlabel="X", ylabel="Y")
	heatmap!(plot_wigner(ρ, 2, 5)..., colormap=:diverging, colorrange=(-.01,.01))
	f
end

plot_wigner(dm(ψ0))

# %%
params = Dict(
	:κ=>.1,
	:η=>1,
	:n=>.1,
	:m=>.1,
	:E=>1,
	:γ=>.1,
	:g=>1
)

# make the timeevolution
times = 0:.1:100
@time tout, ρ = timeevolution.master(times, ψ0, H(;params...), [J(;params...)]);
# add the measurement
# @time tout, ρs = stochastic.master(times, ψ0, H(;params...), [J(;params...)], [C(;params...)], dt=.0001);
#%%
# look at the time evolution
f = Figure(size=fullsize)
a = Axis(f[1,1], ylabel="Q", xlabel="T")
# for (p, l) in zip([ρs, ρ], ["stochastic", "master"])
for (p, l) in zip([ρ], ["master", "stochastic"])
	lines!(times, abs.(expect(Q, p)), label=l)
end
axislegend(position=:lt)
save("../figures/01 time evolution.pdf", f)
f

# %%

f = Figure(size=fullsize)
a = Axis(f[1,1], ylabel="S_Y", xlabel="f")
for (p, l) in zip([ρ, ρs], ["master", "stochastic"])
	freq, S = freq_fft(expect(Y, p), 1/times[2]) 
	lines!(freq/2, abs.(S), label=l)
end
axislegend()
ylims!(a, low=0)
xlims!(a, 0, 2)
save("../figures/01 time evolution spectrum.pdf", f)
f
