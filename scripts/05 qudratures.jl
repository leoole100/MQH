using CairoMakie, QuantumOptics, FFTW, ThreadTools
cd(@__DIR__)
include("00_functions.jl")

# The unites hamiltonian is given:
X = δa + dagger(δa)
Y = im*(dagger(δa)-δa)
Q = δb + dagger(δb)
P = im*(dagger(δb)-δb)

n = dagger(δa)*δa
m = dagger(δb)*δb

H(;Δ=2π, ω=2π, G=0.1) = Δ/2*(X^2 + Y^2) + ω*(Q^2 + P^2) + 2*G*X*Q

J(;κ=1, γ=1,	n=0.1) = √κ * δa + √(γ*(n+1)) * δb + √(γ*n) * dagger(δb)

C(;κ=1, η=1, X=X) = η*√κ*X

# ψ0 = tensor(fockstate(optical_space, 0), fockstate(mechanical_space, 0))
ψ0 = tensor(coherentstate(optical_space, 0), coherentstate(mechanical_space, 0))
ρ0 = dm(ψ0)

function freq_fft(s, fs=1)
	S = fft(s)
	freq = fftfreq(length(s), fs)
	S = S[freq.>0]
	freq = freq[freq.>0]
	return freq, S
end

# %%
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

plot_wigner(ρ0)

# %%
κ=.5
η=.1

# make the timeevolution
times = 0:.1:10
@time tout, ρ = timeevolution.master(times, ψ0, H(), [J(κ=κ)]);
# add the measurement
@time tout, ρs = stochastic.master(times, ψ0, H(), [J(κ=κ)], [C(κ=κ, η=η)], dt=.0001, noise=W);

#%%
# look at the time evolution
f = Figure(size=fullsize)
a = Axis(f[1,1], ylabel="X", xlabel="T")
for (p, l) in zip([ρ, ρs], ["master", "stochastic"])
	lines!(times, abs.(expect(X, p)), label=l)
end
axislegend()
save("../figures/05 time evolution.pdf", f)
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
save("../figures/05 frequency.pdf", f)
f

# %%
# let's sweep the power
function simulation(;κ=.5, η=.1, G=1)
	times = 0:.1:10
	tout, ρ = timeevolution.master(times, ψ0, H(G=G), [J(κ=κ)]);
	tout, ρs = stochastic.master(times, ψ0, H(G=G), [J(κ=κ)], [C(κ=κ, η=η)], dt=.0001);
	return [ρ, ρs]
end

drives = logrange(.5, 10, length=10)
# @time sims = hcat([simulation(E=E) for E in drives]...);
@time sims = hcat(tmap(G->simulation(G=G), drives)...);


# %%
# calculate the peak intensity
function S_max(s, O=X)
	freq, S = freq_fft(expect(O, s))
	return maximum(abs.(S))
end

S = S_max.(sims)

f = Figure(size=fullsize)
a = Axis(f[1,1], ylabel="Sₓ max", xlabel="drive strength E", 
	# xscale=log10, yscale=log10
)
scatterlines!(drives, S[1,:], label="master")
scatterlines!(drives, S[2,:], label="stochastic")
axislegend(position = :lt)
save("../figures/05 power.pdf", f)
f
