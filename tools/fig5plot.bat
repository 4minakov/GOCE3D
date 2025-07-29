gmt begin ../fig/fig5 pdf
  gmt subplot begin 4x2 -Fs8c -M0 -A+JTL -R-56/14/57/80 -Js-20/90/12c/50 -Ba
    gmt subplot set
	gmt makecpt -Croma -Ic -T-1/1/.2 -Z 
    gmt grdimage ../data/Trr_syn.nc -E150 -V -B+tTrr
	rem gmt psscale -B.5x -By+lE
    gmt subplot set
	gmt makecpt -Croma -Ic -T-1/1/.2 -Z 
    gmt grdimage ../data/Trr_rec_syn.nc -E150 -V
	rem gmt psscale -B.5x -By+lEp
    gmt subplot set
	gmt makecpt -Croma -Ic -T-1/1/.2 -Z 
    gmt grdimage ../data/Tzz_syn.nc -E150  -V -B+tTzz
	rem gmt psscale -B.5x -By+lE
    gmt subplot set
	gmt makecpt -Croma -Ic -T-1/1/.2 -Z 
    gmt grdimage ../data/Tzz_rec_syn.nc -E150 -V
	rem gmt psscale -B.5x -By+lE
	gmt subplot set
	gmt makecpt -Croma -Ic -T-1/1/.2 -Z 
    gmt grdimage ../data/Txx_syn.nc -E150 -V -B+tTxx
	rem gmt psscale -B.5x -By+lE
	gmt subplot set
	gmt makecpt -Croma -Ic -T-1/1/.2 -Z 
    gmt grdimage ../data/Txx_rec_syn.nc -E150 -V
	rem gmt psscale -B.5x -By+lE
	gmt subplot set
	gmt makecpt -Croma -Ic -T-1/1/.2 -Z 
    gmt grdimage ../data/Tzx_syn.nc -E150 -V -B+tTzx
	gmt psscale -B.5x -By+lE
	gmt subplot set
	gmt makecpt -Croma -Ic -T-1/1/.2 -Z 
    gmt grdimage ../data/Tzx_rec_syn.nc -E150 -V
	gmt psscale -B.5x -By+lE
  gmt subplot end
gmt end show
