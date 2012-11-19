#!/usr/bin/gnuplot -persist
reset
set terminal postscript eps enhanced color solid
set encoding koi8r
set xlabel "Number of mpi processes"
set ylabel "Mean time, s "
set ytics nomirror
set autoscale
set output "mean_time.ps"

#set title ""
set key reverse Left outside
set grid
set size ratio 0.75
set style line 1 lt 1 pt 5
set style line 2 lt 3 pt 7

plot "performance_data.txt" u 2:3:(column(-2)) t "" w linespoints linestyle 2 lc variable axis x1y2   
