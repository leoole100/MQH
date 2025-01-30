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
fullsize = (4inch, (4/1.33)inch)
halfsize = (2inch, (4/1.33)inch)


using QuantumOptics, LinearAlgebra, StatsBase

function measure(C::Operator, ψ::Ket)
	eigenvals, eigenvecs = eigen(dense(C).data)
	
	probs = [abs2(dagger(ψ)*Ket(basis,eigenvecs[:, i])) for i in 1:length(eigenvals)]
	
	outcome_index = sample( 1:length(eigenvals), Weights(probs))
	outcome_value = eigenvals[outcome_index]
	ψ_post = normalize(Ket(basis, eigenvecs[:, outcome_index]))
	
	return outcome_value, ψ_post
end

function expect_uncert(O::Operator, ψs)
	expectation = expect(O, ψs)
	uncertainty = sqrt.(expect(O*O, ψs) .- expectation.^2)
	return expectation, uncertainty
end

function plot_sim(measured, color; O=sigmaz(basis), uncertainty=true)
	obs = real.(expect_uncert(O, measured[2]))
	if uncertainty
		band!(measured[1], obs[1]-obs[2], obs[1]+obs[2], alpha=.2, color=color)
	end
	return lines!(measured[1], obs[1], color=color)
end
