
set xlabel 'h' ; set ylabel 'ninf'
set logscal xy
set format x "10^{%L}"
set format y "10^{%L}"
set grid xtics
set grid mxtics 
set grid ytics
set grid mytics 
set autoscale fix
set pointsize 3
set title 'Norme infini en fonction du pas h log/log'
set nokey
plot 'ordre.out' using 3:2 with points lc rgb "red"
set terminal png
set output "ordre_h_ninf.png"
replot


