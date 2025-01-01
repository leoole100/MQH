# make figure style
using Makie
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
halfsize = (2inch, (2/1.618)inch)
fullsize = (4inch, (4/1.618)inch)