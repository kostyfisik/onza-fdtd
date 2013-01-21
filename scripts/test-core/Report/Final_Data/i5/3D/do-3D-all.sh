#!/usr/bin/gnuplot -persist
reset
set terminal postscript enhanced color
set encoding koi8r
set xlabel "Task spatial size" font "Times-Roman,35"
set ylabel "Evaluation time per cell, s" font "Times-Roman,35"
#set autoscale
set xrange[20:140]
set yrange[0:1e-7]
set key horiz
set key box
set key ins center bottom
set grid
set size ratio 0.75
set style line 1 lt 1 pt 5 lw 3 linecolor rgb "red"
set style line 2 lt 1 pt 7 lw 3 linecolor rgb "blue"
set style line 3 lt 1 pt 9 ps 1.0 lw 3 linecolor rgb "green"
set style line 4 lt 1 pt 13 lw 3 linecolor rgb "black"
set output '| ps2pdf - 3D_i5.pdf'

plot "3D_test_performance.txt" u 1:2 t "Stencil" w linespoints linestyle 1 axis x1y1, \
"3D_test_performance.txt" u 1:4 t "Range" w linespoints linestyle 2 axis x1y1, \
"3D_test_performance.txt" u 1:6 t "C" w linespoints linestyle 3 axis x1y1, \
"3D_test_performance.txt" u 1:8 t "Fortran" w linespoints linestyle 4 axis x1y1
