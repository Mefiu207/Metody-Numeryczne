#set term postscript size 20cm, 40cm color enhanced solid font "Helvetica, 10"

#
#
#   plik z danymi: x, y, u(x,y)
#
#
reset
set term png size 2600, 3000 enhanced font "Helvetica, 24" lw 2
#set size 1,1
#set size square

set lmargin  7.5
set bmargin  1.5
set rmargin  9.5
set tmargin  1.0


set out "u_t_D0_D01.png"
set grid
set multiplot layout 5,2

set grid
#set object 1 rect from 0,0 to 0.5,0.55 fs solid noborder fc rgb "black" front # przeszkoda

#------------------------------------------------------------------------

set pm3d
set view map
unset surface

#set palette rgbformulae 7,5,15 # paleta kolorow
set palette rgbformulae 22,13,-31
#set palette rgbformulae 33,13,10

#set palette defined (-2 "blue", 0 "white", 16 "red")


#------------------------------------------------------------------------
set cbrange[0:18]
set title 'D=0 u(x,y) it=2000, t=0.1139'
set xlabel 'x'
set ylabel 'y'

splot [0:400][0:90] "u_t_d0.dat" u 1:2:3 i 1 notitle


#------------------------------------------------------------------------
set cbrange[0:4]
set title 'D=0.1 u(x,y) it=2000, t=0.1139'
set xlabel 'x'
set ylabel 'y'

splot [0:400][0:90] "u_t_d01.dat" u 1:2:3 i 1 notitle

#------------------------------------------------------------------------

#set title 'u(x,y) it=3000, t='
#set xlabel 'x'
#set ylabel 'y'

#splot [0:400][0:90] "u_t.dat" u 1:2:3 i 2 notitle

#------------------------------------------------------------------------

#set title 'u(x,y) it=, t='
#set xlabel 'x'
#set ylabel 'y'

#splot [0:400][0:90] "u_t.dat" u 1:2:3 i 3 notitle

#------------------------------------------------------------------------
set cbrange[0:18]
set title 'D=0 u(x,y) it=4000, t=0.2279'
set xlabel 'x'
set ylabel 'y'

splot [0:400][0:90] "u_t_d0.dat" u 1:2:3 i 3 notitle

#------------------------------------------------------------------------
set cbrange[0:2.5]
set title 'D=0.1 u(x,y) it=4000, t=0.2279'
set xlabel 'x'
set ylabel 'y'

splot [0:400][0:90] "u_t_d01.dat" u 1:2:3 i 3 notitle

#------------------------------------------------------------------------

#set title 'u(x,y) it=5000, t='
#set xlabel 'x'
#set ylabel 'y'

#splot [0:400][0:90] "u_t.dat" u 1:2:3 i 5 notitle

#------------------------------------------------------------------------
set cbrange[0:18]
set title 'D=0 u(x,y) it=6000, t=0.3418'
set xlabel 'x'
set ylabel 'y'

splot [0:400][0:90] "u_t_d0.dat" u 1:2:3 i 5 notitle

#------------------------------------------------------------------------
set cbrange[0:1.2]
set title 'D=0.1 u(x,y) it=6000, t=0.3418'
set xlabel 'x'
set ylabel 'y'

splot [0:400][0:90] "u_t_d01.dat" u 1:2:3 i 5 notitle

#------------------------------------------------------------------------

#set title 'u(x,y) it=7000, t='
#set xlabel 'x'
#set ylabel 'y'

#splot [0:400][0:90] "u_t.dat" u 1:2:3 i 7 notitle

#------------------------------------------------------------------------
set cbrange[0:18]
set title 'D=0 u(x,y) it=8000, t=0.4558'
set xlabel 'x'
set ylabel 'y'

splot [0:400][0:90] "u_t_d0.dat" u 1:2:3 i 7 notitle

#------------------------------------------------------------------------
set cbrange[0:0.8]
set title 'D=0.1 u(x,y) it=8000, t=0.4558'
set xlabel 'x'
set ylabel 'y'

splot [0:400][0:90] "u_t_d01.dat" u 1:2:3 i 7 notitle
#------------------------------------------------------------------------

#set title 'u(x,y) it=9000, t='
#set xlabel 'x'
#set ylabel 'y'

#splot [0:400][0:90] "u_t.dat" u 1:2:3 i 8 notitle

#------------------------------------------------------------------------
set cbrange[0:18]
set title 'D=0 u(x,y) it=10000, t=0.5697'
set xlabel 'x'
set ylabel 'y'

splot [0:400][0:90] "u_t_d0.dat" u 1:2:3 i 9 notitle

#------------------------------------------------------------------------
set cbrange[0:0.5]
set title 'D=0.1 u(x,y) it=10000, t=0.5697'
set xlabel 'x'
set ylabel 'y'

splot [0:400][0:90] "u_t_d01.dat" u 1:2:3 i 9 notitle

#------------------------------------------------------------------------



unset multiplot
reset

