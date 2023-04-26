set terminal postscript eps enhanced color font 'Helvetica,10' size 7.5cm,7.5cm
set output 'land_survey.eps'

#set terminal pngcairo size 800,800
# set xtics font "Helvetica,15"
# set ytics font "Helvetica,15"
# set key font "Helvetica,15"
# set title font 'Helvetica,20'
#set output '| display png:-'


#set multiplot layout 2,1
set grid back
# set xrange [-3000:3000]
# set yrange [-3000:3000]
set title "Acquisition geometry"
set xlabel 'X(m)'
set ylabel 'Y(m)'

set label '1' at -1600.00000,-1600.00000 font ',16'
set label '2' at -533.333374,-1600.00000 font ',16'
set label '3' at  533.333252,-1600.00000 font ',16'
set label '4' at  1600.00000,-1600.00000 font ',16'
set label '5' at -1600.00000,-533.333374 font ',16'
set label '6' at -533.333374,-533.333374 font ',16'
set label '7' at 533.333252,-533.333374 font ',16'
set label '8' at  1600.00000,-533.333374 font ',16'
set label '9' at  -1600.00000,533.333252 font ',16'
set label '10' at  -533.333374,533.333252 font ',16'
set label '11' at  533.333252,533.333252 font ',16'
set label '12' at   1600.00000,533.333252 font ',16'
set label '13' at  -1600.00000,1600.00000 font ',16'
set label '14' at  -533.333374,1600.00000 font ',16'
set label '15' at  533.333252,1600.00000 font ',16'
set label '16' at   1600.00000,1600.00000 font ',16'


plot "receivers.txt" using 1:2 with points ps 1 pt 8 title 'Receivers', \
"sources.txt" using 1:2 with points ps 1 pt 7 title 'Transmitters' 


# plot "receivers.txt" using 1:2 with lp title 'Receivers', \
# "sources.txt" using 1:2 with points ps 3 pt 7 title 'Transmitters' 

