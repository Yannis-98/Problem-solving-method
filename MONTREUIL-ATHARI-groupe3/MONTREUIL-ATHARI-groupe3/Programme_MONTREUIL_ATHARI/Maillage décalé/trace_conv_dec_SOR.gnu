
set xlabel 'x' ; set ylabel 'U(x)'
set grid xtics
set grid mxtics 
set grid ytics
set grid mytics 
set autoscale fix
set pointsize 1.4
set title 'Solutions U decallé avec SOR en fonction de x'
set nokey
set xrange[-0.25:3.5]
set yrange[-0.2:1.2]
plot 'Solutions_SOR_dec.out' with points lc rgb "blue"
set terminal png
set output "Solutions_SOR_dec.png"
replot


