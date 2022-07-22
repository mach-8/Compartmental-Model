#!/usr/local/bin/gnuplot

#########################################################################################
#
#   BUBBLE SIZE DISTRIBUTION
#
#########################################################################################

reset

set terminal pdfcairo enhanced color dashed rounded size 16 cm, 9.6 cm font "Alegreya, 14"

set output 'size_distribution.pdf'

load 'color.pal'

load 'xyborder.cfg'

load 'grid.cfg'

filenames(n) = sprintf('../OUTPUTS/QGT-%d.dat',n)
legend(val) = sprintf('Q_{gt} = %.2f m^3/s',val)
linestyle(n) = n/10

set xlabel 'volume equivalent bubble diameter, d_{bi} (mm)'
set ylabel 'probability density'
set key right top
plot for [i = 10:50:10] filenames(i) every ::1 u 1:3 w l ls linestyle(i) t legend(i/100.0)

#########################################################################################
#
#   BUBBLE SLIP VELOCITY DISTRIBUTION
#
#########################################################################################
reset

set terminal pdfcairo enhanced color dashed rounded size 16 cm, 9.6 cm font "Alegreya, 14"

set output 'slip_velocity_distribution.pdf'

load 'color.pal'

load 'xyborder.cfg'

load 'grid.cfg'

filenames(n) = sprintf('../OUTPUTS/QGT-%d.dat',n)
legend(val) = sprintf('Q_{gt} = %.2f m^3/s',val)
linestyle(n) = n/10

set xlabel 'bubble slip velocity, u_{si} (m/s)'
set ylabel 'probability density'
set key right top
plot for [i = 10:50:10] filenames(i) every ::1 u 2:3 w l ls linestyle(i) t legend(i/100.0)

#########################################################################################
#
#   HOLDUPS
#
#########################################################################################

reset

set terminal pdfcairo enhanced color dashed rounded size 16 cm, 20 cm font "Alegreya, 14"

set output 'holdups.pdf'

load 'color.pal'

load 'xyborder.cfg'

load 'grid.cfg'

filename = sprintf('../OUTPUTS/summary_results.dat')

set multiplot layout 2,1
    set tmargin 2
    set bmargin 2
    set xlabel 'fresh treat gas velocity, U_{gt} (m/s)'
    set ylabel 'gas holdup (-)'
    set key right bottom
    plot filename every ::1 u 1:6 w l ls 1 t 'bed region', \
    filename every ::1 u 1:8 w l ls 2 t 'freeboard region'

    set tmargin 2
    set bmargin 2
    set ylabel 'bed liquid holdup (-)'
    plot filename every ::1 u 1:7 w l ls 3 t ''
unset multiplot

#########################################################################################
#
#   RECYCLE PAN
#
#########################################################################################

reset

set terminal pdfcairo enhanced color dashed rounded size 16 cm, 20 cm font "Alegreya, 14"

set output 'recycle_pan.pdf'

load 'color.pal'

load 'xyborder.cfg'

load 'grid.cfg'

filename = sprintf('../OUTPUTS/summary_results.dat')

set multiplot layout 2,1
    set tmargin 2
    set bmargin 2
    set xlabel 'fresh treat gas velocity, U_{gt} (m/s)'
    set ylabel ''
    set key right top
    plot filename every ::1 u 1:10 w l ls 1 t 'liquid recycle ratio', \
    filename every ::1 u 1:12 w l ls 2 t 'recycle gas fraction'

    set tmargin 2
    set bmargin 2
    set ylabel 'gas-liquid separation efficiency'
    plot filename every ::1 u 1:9 w l ls 3 t ''
unset multiplot