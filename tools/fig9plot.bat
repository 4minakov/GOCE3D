gmt begin ../fig/fig9new pdf
  gmt subplot begin 2x2 -Fs8c -M0 -A+JTL -R-60/20/55/82 -Js-20/90/12c/50 -Ba
    gmt subplot set
	rem gmt grd2cpt Trr_topofull.nc -Croma -L-3/3
	gmt makecpt -Croma -T-3/3/.5
    gmt grdimage ../data/Trr_topofull.nc -E150 -V
	gmt psscale -B1x -By+lE
	gmt coast -Di -W0.5p -A1000 -N1/0.5p
    gmt subplot set
	rem gmt grd2cpt Trr_sed.nc -Croma -L-1/1
	gmt makecpt -Croma -T-3/1/.5
    gmt grdimage ../data/Trr_sed.nc -E150 -V
	gmt psscale -B.5x -By+lE
	gmt coast -Di -W0.5p -A1000 -N1/0.5p
    gmt subplot set
	rem gmt grd2cpt Trr_moho.nc -Croma -L-5/5
	gmt makecpt -Croma -T-5/5/1
    gmt grdimage ../data/Trr_moho.nc -E150  -V
	gmt psscale -B1x -By+lE
	gmt coast -Di -W0.5p -A1000 -N1/0.5p
    gmt subplot set
	rem gmt grd2cpt Trr_therm.nc -Croma -L-3/1
	gmt makecpt -Croma -T-3/0.5/.5
    gmt grdimage ../data/Trr_therm.nc -E150 -V
	gmt psscale -B.5x -By+lE
	gmt coast -Di -W0.5p -A1000 -N1/0.5p
  gmt subplot end
gmt end show
