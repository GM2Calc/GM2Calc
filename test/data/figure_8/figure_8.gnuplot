set terminal pdf color enhanced
set output 'figure_8.pdf'
set logscale y 10
set grid
set key box center right
data = 'figure_8.txt'

set style line 1 lc '#0000ff'
set style line 2 lc '#ff0000'
set style line 3 lc '#00aa00'

set ylabel '|a_μ / ζ_fζ_l|'
set xlabel 'M_H / GeV'

plot [0:500] [5e-15:1e-10] \
     data u 1:(abs($2)) w lines t 'u' ls 1, \
     data u 1:(abs($3)) w lines t 'd' ls 2, \
     data u 1:(abs($4)) w lines t 'l' ls 3

set xlabel 'M_A / GeV'

plot [0:500] [5e-15:1e-10] \
     data u 1:(abs($5)) w lines t 'u' ls 1, \
     data u 1:(abs($6)) w lines t 'd' ls 2, \
     data u 1:(abs($7)) w lines t 'l' ls 3

set xlabel 'M_{H^±} / GeV'

plot [0:500] [1e-16:1e-10] \
     data u 1:(abs($8))  w lines t 'u' ls 1, \
     data u 1:(abs($9))  w lines t 'd' ls 2, \
     data u 1:(abs($10)) w lines t 'l' ls 3

set ylabel '|a_μ / ηζ_l|'
set xlabel 'M_H / GeV'

plot [0:500] [1e-14:1e-10] \
     data u 1:(abs($11)) w lines t 'u' ls 1, \
     data u 1:(abs($12)) w lines t 'd' ls 2, \
     data u 1:(abs($13)) w lines t 'l' ls 3
