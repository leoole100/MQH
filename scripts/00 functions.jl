cd(@__DIR__)

# make figure style
using CairoMakie
inch = 96
pt = 4/3
cm = inch/2.54
update_theme!(
	fontsize=10pt,
	px_per_unit = 10,
	fonts = (;
		regular="Roboto Regular",
		bold="Roboto Bold",
	),
)
halfsize = (2inch, (2/1.33)inch)
fullsize = (4inch, (4/1.33)inch)


# Define the basis
using QuantumOptics

optical_space = FockBasis(10)
mechanical_space = FockBasis(5)
hilbert_space = optical_space ⊗ mechanical_space

# and the commonly used operators
δa = destroy(optical_space) ⊗ one(mechanical_space)
δb = one(optical_space) ⊗ destroy(mechanical_space)

# The unites hamiltonian is given:
X = δa + dagger(δa)
Y = im*(dagger(δa)-δa)
Q = δb + dagger(δb)
P = im*(dagger(δb)-δb)

n = dagger(δa)*δa
m = dagger(δb)*δb

# H(;ω=2π, Δ=-ω, G=0.1, kwargs...) = Δ/2*(X^2 + Y^2) + ω*(Q^2 + P^2) + 2*G*X*Q
H(;ω=2π, Δ=-ω, g=0.1, E=1, kwargs...) = Δ*n+ω*m-g*n*Q+E*X
# J(;κ=1, γ=1,	m=0.1, kwargs...) = √κ * δa + √(γ*(m+1)) * δb + √(γ*m) * dagger(δb)
J(;κ=1, γ=1,	n=0, m=0.1, kwargs...) = √(κ*(n+1)) * δa + √(κ*n) * dagger(δa) + √(γ*(m+1)) * δb + √(γ*m) * dagger(δb)
# C(;κ=1, η=1, X=X, kwargs...) = η*√κ*X
C(;κ=1, η=1, X=X, kwargs...) = η*√κ*δa

# ψ0 = tensor(fockstate(optical_space, 0), fockstate(mechanical_space, 0))
ψ0 = tensor(coherentstate(optical_space, 0), coherentstate(mechanical_space, 0))

using FFTW
function freq_fft(s, fs=1)
	S = fft(s)
	freq = fftfreq(length(s), fs)
	S = S[freq.>0]
	freq = freq[freq.>0]
	return freq, S
end