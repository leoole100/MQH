include("00 functions.jl")

# %%
params = Dict(
	:κ=>0.1,
	:η=>0.1,
	:n=>0.1,
	:m=>0.1,
	:γ=>0.1,
	:g=>0.1,
)
drives = logrange(1e-3, 1e2, length=10)
C(;κ=1, η=1, X=X, kwargs...) = η*√κ*Q

# let's sweep the power
function simulation(;E=1., params=params)
	params = merge(params, Dict(:E=>E))
	times = 0:.1:100
	tout, ρ = timeevolution.master(times, ψ0, H(;params...), [J(;params...)]);
	ρs = ρ
	tout, ρs = stochastic.master(times, ψ0, H(;params...), [J(;params...)], [C(;params...)], dt=.0001);
	return [ρ, ρs]
end

# @time sims = hcat([simulation(E=E) for E in drives]...);
using ThreadTools
@time sims = hcat(tmap(E->simulation(;E=E), drives)...); # takes 10 min


# %%
# calculate the peak intensity
function S_max(s, O=C(;params...))
	freq, S = freq_fft(expect(O, s))
	return maximum(abs.(S))
end

S = S_max.(sims) 

f = Figure(size=fullsize)
a = Axis(f[1,1], ylabel="S_C max", xlabel="drive strength E, α", 
	xscale=log10, yscale=log10
)
scatterlines!(drives, S[1,:], label="master")
scatterlines!(drives, S[2,:], label="stochastic")
function analytical(α; params=params, ω=1)
	return params[:κ]*params[:γ]*abs(1-2im*ω/params[:κ])/(8*params[:η]*params[:γ]*4*params[:g]^2*α^2) +
	2*params[:γ]*(2π)^2/abs((2π)^2-ω^2+1im*ω*params[:γ])^2 *
	4*params[:g]^2*α^2/(params[:κ]*params[:γ]*(1-2im*ω/params[:κ])^2)
	
end
d = logrange(minimum(drives), maximum(drives), length=100)
lines!(d,	abs.(analytical.(d)), label="analytical solution", color=:gray)

axislegend(position = :ct, framevisible=false)
save("../figures/02 power.pdf", f)
f
