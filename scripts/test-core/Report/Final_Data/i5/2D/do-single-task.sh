#!/usr/bin/gnuplot -persist
reset
set terminal postscript enhanced color
set encoding koi8r
set xlabel "Blocks' size" font "Times-Roman,35"
set ylabel "Evaluation time per cell, s" font "Times-Roman,35"
#set autoscale
set xrange[0:520]
set yrange[1e-8:7e-8]
set key horiz
set key box
set key ins center top
set grid
set size ratio 0.75
set style line 1 lt 1 pt 5 ps 1.0 lw 3 linecolor rgb "red"
set style line 2 lt 1 pt 7 lw 3 linecolor rgb "blue"
set style line 3 lt 1 pt 0 lw 3 linecolor rgb "red"
set style line 4 lt 1 pt 0 lw 3 linecolor rgb "blue"
set output '| ps2pdf - 1024_blocks_i5.pdf'

plot "results_size_1024x1024.plot" u 1:7 t "Stencil" w linespoints linestyle 3 axis x1y1, \
"results_size_1024x1024.plot" u 1:8 t "Stencil (blocks)" w linespoints linestyle 1 axis x1y1, \
"results_size_1024x1024.plot" u 1:9 t "Range" w linespoints linestyle 4 axis x1y1, \
"results_size_1024x1024.plot" u 1:10 t "Range (blocks)" w linespoints linestyle 2 axis x1y1
