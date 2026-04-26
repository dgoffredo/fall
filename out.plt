set title 'Displacement versus Time as You Fall Through the Center of the Earth'
set ylabel 'Distance from Center of the Earth (kilometers)'
set xlabel 'Time Elapsed (minutes)'

set mxtics
set mytics
set grid

plot 'out' using ($1/60):($2/1000) notitle
