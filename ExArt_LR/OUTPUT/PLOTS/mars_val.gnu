set term pdf col enh dashed size 17cm,17cm

set out "validation_mars.pdf"

set border lw 3
set origin 0,0
set size 1,1.5
set key box

set termoption dashed

set log y

set multiplot layout 2,2

set grid

#################################                             |
#################################                            c|

set origin 0,0
set size 0.5,0.5
set key left top

set xlabel "Flux [W/m^2]"
set ylabel "Pressure [mbar]"

a="Our Code"
b="AMES 1D/DISORT"
c="AMES 1D/AMES
d="KDM/DISORT"
e="KDM/KDM"
f="AMES 1D/CLOUDS"

set ytics (0.01,0.1,1.0,10.0,100.0,500.0)
set xtics (160,170,180,190,200,210,220)
set label "DOWN SOLAR" at 170,0.3 font ",18" 
p [160:230][600:0.01] 'SOL_ExArt_mars.out' u 2:($1/100) w l lc 1 lw 4 t a,\
               'SOL_AMES_mars.out' u 2:($1/100) w l lc 4 lt 4 lw 4 t b,\
               'SOL_AMES_OR_mars.out' u 2:1 w l lc 2 lw 4 t c,\
               'SOL_AMES_OR_mars_CL.out' u 2:1 w l lc 6 lw 4 t f,\
               'SOL_MIS_mars.out' u 2:($1/100) w l lc 3 lt 2 lw 4 t d,\
               'SOL_MIS_OR.out' u 2:1 w l lc 5 lt 2 lw 5 t e

#################################                            c|
#################################                             |

set size 0.5,0.5
set origin 0.,0.5
set xlabel ""

set ytics (0.01,0.1,1.0,10.0,100.0,500.0)
set xtics (20,25,30,35,40,45,50,55,60,65,70)
unset label
set label "UP SOLAR" at 5,0.3 font ",18" 

p [:70][510:0.01] 'SOL_ExArt_mars.out' u 4:($1/100) w l lc 1 lw 4 t a,\
               'SOL_AMES_mars.out' u 4:($1/100) w l lw 4 lc 4 lt 4 t b,\
               'SOL_AMES_OR_mars.out' u 3:1 w l lc 2 lw 4 t c,\
               'SOL_AMES_OR_mars_CL.out' u 3:1 w l lc 6 lw 4 t f,\
               'SOL_MIS_mars.out' u 4:($1/100) w l lc 3 lt 2 lw 4 t d,\
               'SOL_MIS_OR.out' u 3:1 w l lc 5 lt 2 lw 5 t e
               
#################################                             |c
#################################                             |
unset ylabel
set size 0.5,0.5
set origin 0.5,0.5
set xlabel ""

unset label
set label "UP IR" at 170,0.3 font ",18" 
set key right top

set xtics (140,160,180,200,220)
p [120:][520:0.01] 'IR_ExArt_mars.out' u 4:($1/100) w l lc 1 lw 4 t a,\
               'IR_AMES_mars.out' u 4:($1/1e2) w l lw 4 lc 4 lt 4 t b,\
               'IR_AMES_OR_mars.out' u 3:1 w l lc 2 lw 4 t c,\
               'IR_AMES_OR_mars_CL.out' u 3:1 w l lc 6 lw 4 t f,\
               'IR_MIS_mars.out' u 4:($1/100) w l lc 3 lw 4 lt 2 t d,\
               'IR_MIS_OR.out' u 3:1 w l lc 5 lw 5 lt 2 t e

#################################                             |c
#################################                             |

set size 0.5,0.5
set origin 0.5,0.
set xlabel "Flux [W/m^2]"


unset label
set label "DOWN IR" at 50,0.3 font ",18" 

set xtics (20,40,60,80,100,120,140)
p [:150][510:0.01] 'IR_ExArt_mars.out' u 2:($1/100) w l lw 4 lc 1 t a,\
               'IR_AMES_mars.out' u 2:($1/1e2) w l lc 4 lt 4 lw 4 t b,\
               'IR_AMES_OR_mars.out' u 2:1 w l lw 4 lc 2 t c,\
               'IR_AMES_OR_mars_CL.out' u 2:1 w l lw 4 lc 6 t f,\
               'IR_MIS_mars.out' u 2:($1/100) w l lw 4 lc 3 lt 2 t d,\
               'IR_MIS_OR.out' u 2:1 w l lw 5 lc 5 lt 2 t e

