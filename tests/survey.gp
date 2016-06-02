# usage:
# gnuplot -e "file='survey.tfs'; title='sequ'" survey.gp

set grid
set xlabel 'X'
set ylabel 'Z'
set size square
set offsets 1, 1, 1, 1
plot file u 6:8 w lp t title

pause -1
