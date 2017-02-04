rm -f  sbend.tfsone
rm -f  table_sbend_ptc.out
rm -f  table_sbend_madng.out
rm -f  table_sbend.out
rm -f  sbend_ptc.out
rm -f  sbend_madng.out
rm -f  sbend.eps

madxl sbend_ptc.madx > sbend_ptc.out
./mad all.mad -v -p TestTrack.testTrackSBEND > sbend_madng.out

grep -i "PTC_SBEND"   sbend_ptc.out   > table_sbend_ptc.out
grep -i "MADNG_SBEND" sbend_madng.out > table_sbend_madng.out
paste table_sbend_ptc.out table_sbend_madng.out > table_sbend.out

gnuplot --persist <<EOF
  set term post eps enh col
  set out 'sbend.eps'
  set title "SBEND PTC vs MAD-NG, Yoshida6, nst=5 "
  set ylabel 'x'
  set y2label 'x(PTC - MAD-NG)'
  set xlabel 'k0, sbend strength'
  set grid x y
  set ytics nomirror
  set y2tics
  set tics out
  set autoscale  y
  set autoscale y2

  plot 'table_sbend.out' u 3:4 t 'PTC' w lp lw 2 axes x1y1,\
  'table_sbend.out' u 12:13 t 'MAD-NG' w lp lw 2 axes x1y1,\
  'table_sbend.out' u 3:(\$4-\$13) t 'X(PTC - MAD-NG)' w lp lw 2 axes x1y2

  set ylabel 'px'
  set y2label 'px(PTC - MAD-NG)'
  set xlabel 'k0, sbend strength'
  set grid x y
  set ytics nomirror
  set y2tics
  set tics out
  set autoscale  y
  set autoscale y2

  plot 'table_sbend.out' u 3:5 t 'PTC' w lp lw 2 axes x1y1,\
  'table_sbend.out' u 12:14 t 'MAD-NG' w lp lw 2 axes x1y1,\
  'table_sbend.out' u 3:(\$5-\$14) t 'PX(PTC - MAD-NG)' w lp lw 2 axes x1y2

  set ylabel 'y'
  set y2label 'y(PTC - MAD-NG)'
  set xlabel 'k0, sbend strength'
  set grid x y
  set ytics nomirror
  set y2tics
  set tics out
  set autoscale  y
  set autoscale y2

  plot 'table_sbend.out' u 3:6 t 'PTC' w lp lw 2 axes x1y1,\
  'table_sbend.out' u 12:15 t 'MAD-NG' w lp lw 2 axes x1y1,\
  'table_sbend.out' u 3:(\$6-\$15) t 'Y(PTC - MAD-NG)' w lp lw 2 axes x1y2

  set ylabel 'py'
  set y2label 'py(PTC - MAD-NG)'
  set xlabel 'k0, sbend strength'
  set grid x y
  set ytics nomirror
  set y2tics
  set tics out
  set autoscale  y
  set autoscale y2

  plot 'table_sbend.out' u 3:7 t 'PTC' w lp lw 2 axes x1y1,\
  'table_sbend.out' u 12:16 t 'MAD-NG' w lp lw 2 axes x1y1,\
  'table_sbend.out' u 3:(\$7-\$16) t 'PY(PTC - MAD-NG)' w lp lw 2 axes x1y2

  set ylabel 't '
  set y2label 't(PTC - MAD-NG)'
  set xlabel 'k0, sbend strength'
  set grid x y
  set ytics nomirror
  set y2tics
  set tics out
  set autoscale  y
  set autoscale y2

  plot 'table_sbend.out' u 3:(\$8*(-1)) t 'PTC' w lp lw 2 axes x1y1,\
  'table_sbend.out' u 12:17 t 'MAD-NG' w lp lw 2 axes x1y1,\
  'table_sbend.out' u 3:(\$8+\$17) t 'T(PTC - MAD-NG)' w lp lw 2 axes x1y2

  set ylabel 'pt '
  set y2label 't(PTC - MAD-NG)'
  set xlabel 'k0, sbend strength'
  set grid x y
  set ytics nomirror
  set y2tics
  set tics out
  set autoscale  y
  set autoscale y2

  plot 'table_sbend.out' u 3:9 t 'PTC' w lp lw 2 axes x1y1,\
  'table_sbend.out' u 12:18 t 'MAD-NG' w lp lw 2 axes x1y1,\
  'table_sbend.out' u 3:(\$9-\$18) t 'T(PTC - MAD-NG)' w lp lw 2 axes x1y2
EOF
