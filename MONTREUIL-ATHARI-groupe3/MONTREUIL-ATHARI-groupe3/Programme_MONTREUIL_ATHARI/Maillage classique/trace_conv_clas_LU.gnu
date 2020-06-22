
set xlabel 'x' ; set ylabel 'U(x)'
set grid xtics
set grid mxtics 
set grid ytics
set grid mytics 
set autoscale fix
set pointsize 1.4
set title 'Solutions U classique avec LU en fonction de x'
set nokey
set xrange[-0.25:3.5]
set yrange[-0.2:1.2]
plot 'Solutions_LU.out' with points lc rgb "red"
set terminal png
set output "Solutions_LU.png"
replot


