gmt begin ../fig/fig18gmt pdf
  gmt subplot begin 3x2 -Fs8c -M0 -A+JTL -R-56/14/57/80 -Js-20/90/12c/50 -Ba
    gmt subplot set
	gmt makecpt -Croma -Ic -T-2/2/.5 -Z 
    gmt grdimage ../data/Tzz_real.nc -E150  -V -B+tTzz
	gmt coast -Di -W0.5p -A1000 -N1/0.5p
	rem gmt psscale -B.5x -By+lE
    gmt subplot set
	gmt makecpt -Croma -Ic -T-2/2/.5 -Z 
    gmt grdimage ../data/Tzz_rec_real.nc -E150 -V
	gmt coast -Di -W0.5p -A1000 -N1/0.5p
	rem gmt psscale -B.5x -By+lE
	gmt subplot set
	gmt makecpt -Croma -Ic -T-2/2/.5 -Z 
    gmt grdimage ../data/Txx_real.nc -E150 -V -B+tTxx
	gmt coast -Di -W0.5p -A1000 -N1/0.5p
	rem gmt psscale -B.5x -By+lE
	gmt subplot set
	gmt makecpt -Croma -Ic -T-2/2/.5 -Z 
    gmt grdimage ../data/Txx_rec_real.nc -E150 -V
	gmt coast -Di -W0.5p -A1000 -N1/0.5p
	rem gmt psscale -B.5x -By+lE
	gmt subplot set
	gmt makecpt -Croma -Ic -T-2/2/.5 -Z 
    gmt grdimage ../data/Tzx_real.nc -E150 -V -B+tTzx
	gmt coast -Di -W0.5p -A1000 -N1/0.5p
	gmt psscale -B.5x -By+lE
	gmt subplot set
	gmt makecpt -Croma -Ic -T-2/2/.5 -Z 
    gmt grdimage ../data/Tzx_rec_real.nc -E150 -V
	gmt coast -Di -W0.5p -A1000 -N1/0.5p
	gmt psscale -B.5x -By+lE
  gmt subplot end
gmt end show
