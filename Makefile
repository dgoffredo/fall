fall: fall.cpp density.h Makefile
	c++ -Wall -Wextra -pedantic --std=c++20 -O3 -o $@ $<

density.h: density.csv Makefile
	echo '#pragma once' >$@
	echo '' >>$@
	echo '#include <cmath>' >>$@
	echo '' >>$@
	echo 'struct EarthRecord {' >>$@
	echo '    long double r; // radius in meters' >>$@
	echo '    long double ρ; // density in kilograms per cubic meter' >>$@
	echo '    long double μ; // Δρ/Δr relative to the previous record' >>$@
	echo '    long double m; // mass contained at this radius, in kilograms' >>$@
	echo '};' >>$@
	echo '' >>$@
	echo 'constexpr EarthRecord earth_table[] = {' >>$@
	awk -F, '(NR>2) { print "    {.r="$$1", .ρ="$$2", .μ="$$5", .m="$$6"}," }' $< | sed 's,#DIV/0!,NAN,g' >>$@
	echo '};' >>$@

out: fall
	./$< | head -10000 >$@

plot: fall out out.plt
	gnuplot out.plt -

