
set xlabel 'h' ; set ylabel 'const'
set logscal x
set autoscale y
set format y "10^{%L}"
set grid xtics
set grid mxtics 
set grid ytics
set grid mytics 
set autoscale fix
set pointsize 3
set title 'Constante en fonction du pas h lin/log'
set nokey
set yrange[0:1]
plot 'ordre.out' using 2:4 with points lc rgb "red"
set terminal png
set output "ordre_h_const.png"
replot


