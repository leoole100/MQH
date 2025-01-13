include("00_functions.jl")

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
using ThreadTools
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
save("../figures/02 power.pdf", f)
f
