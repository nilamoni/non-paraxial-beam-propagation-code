#	axes
	set xlabel "X/A"
	set ylabel "Z/L"

#	axis ranges
	r_min			=	-2 	
	r_max			=	2
	r_range			=	abs(r_min - r_max)
	dr				=	1e-2/4
	lr				=	r_range/dr
	lr2				=	lr/2

	
	z_max 			=	2.0
	dz				=	1e-4
	lz				=	z_max/dz	

	set yrange[0:z_max]
	set xrange[r_min:r_max]



set terminal png size 1024,1024 enhanced font "Helvetica,20"
set output '0_8w_f_X_0_0_2.png'
set size  square
unset key

set ytics out 0.5
set xtics out 0.5

set title "W_0 = 0.8W_f, X_0/A = 0.2"


skipr = 1
skipz = int(lz/1000)


plot "longitudinal_profille.dat" u (($1-lr2)*dr*skipr):(($2 - 1)*dz*skipz):3  matrix w image,\
0.5 w l dt 2 lc "yellow",\
1.0 w l dt 2 lc "yellow"
#0.25 w l dt 2 lc "yellow",\
#0.75 w l dt 2 lc "yellow",\
1.0 w l dt 2 lc "yellow"


2