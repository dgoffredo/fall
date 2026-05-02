fall: fall.cpp density.h
	c++ -Wall -Wextra -pedantic --std=c++20 -g -O3 -o $@ $<

jsontime: jsontime.cpp
	c++ -Wall -Wextra -pedantic --std=c++20 -O2 -o $@ $<

density.h: density.csv
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

plot: out.plt fall out
	gnuplot $< -

out.png: out.plt fall out
	gnuplot \
		-e 'set terminal png truecolor nocrop enhanced butt size 1920,1080 font "arial,18.0"' \
		-e 'set output "$@"' \
		$<

.PHONY: clean
clean:
	rm -f out.png plot out density.h fall
