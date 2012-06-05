set terminal png nocrop size 1270,600
set output "calc1.png"
set xlabel "X coord"
set ylabel "Ez"
#set yrange [-0.1:0.1]
set yrange [-1:1]
#set linestyle 1 lt 2 lw 3
#set key below 
set title "calc1"
#
plot	"calc1" using 1:2 title 'Ez(x)' w l
#,\
#	"calc1" using 1:($4/$3) title 'transmission' w l
#	"calc1" using 1:4 title 'pmlRight' w l,\
#	"calc1" using 1:5 title 'pmlTop' w l
#	"calc1" using 1:($4/$3) title 'WG-diff' w l
#	"calc1" using 1:26 title 'ver/bl'w l lw 3,\
#	"calc1" using 1:25 title 'hor/bl'w l

#8 tr/bl	fbl	fb	fbr	cr	ftr	ft	ftl	cl	fbl/bl	fb/bl	fbr/bl	cr/bl	#ftr/bl	ft/bl	 ftl/bl	cl/bl	25sum_hor_norm	sum_ver_norm
