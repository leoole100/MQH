include("00 functions.jl")
using QuantumOptics


space = FockBasis(20)

a = destroy(space)
at = dagger(a)

x = a + at
p = im*(at-a)
H() = 1/2*p^2 + 1/2*x^2
ψ0 = coherentstate(space, 1)

dt = 0.001
η = .8
O = x

# C(η=1) = η*x 
C(η=1) = η*x

# time evolution
times = 0:0.1:10

tout, free = timeevolution.schroedinger(times, ψ0, H())
@time measured = [stochastic.master(times, ψ0, H(), [], [C(η)], dt=dt) for i in 1:10];

#%%

f = Figure(size=fullsize)
a = Axis(f[1,1], xlabel="t", ylabel="Re <x>", title="dt = $(dt), η = $(η)")
lines!(times, real.(expect(O, free)), label="free", color=:black)

for (t, m) in measured
	if length(t)==length(times)
		lines!(t, real.(expect(O, m)), color=:green, alpha=.25, label="measured")	
	else
		lines!(t, real.(expect(O, m)), color=:red, alpha=.25, label="instability detected")	
	end
end

axislegend(padding=0, merge=true, position=:lb, framevisible=false)
ylims!(-10, 10)
save("../figures/03 harmonic oscillator.pdf", f)
f
