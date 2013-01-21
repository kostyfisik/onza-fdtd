set terminal png nocrop size 1270,600
set ylabel "Processing Time, s"
set xlabel "Cache line size, elements"
set autoscale
set output "time_name.png"

set grid
set size ratio 0.75
set style line 1 lt 1 pt 5
set style line 2 lt 3 pt 7
set label "name" at graph 0.3, graph 0.9

plot "name" u 1:2 t "StdC" w lines, \
"name" u 1:3 t "StdC with blocks" w linespoints linestyle 1 axis x1y1 
  
