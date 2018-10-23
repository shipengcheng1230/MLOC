#!/bin/bash
gmt gmtset PS_MEDIA letter
gmt gmtset PROJ_LENGTH_UNIT cm
gmt gmtset PS_PAGE_ORIENTATION portrait
psfile=jsa4.1_example_tt7_local_regional_s.ps
projection=-JX18/12
region=-R0/15/-120/20
gmt psbasemap $region $projection -Bxa5f1+l'Epicentral Distance (deg)' -Bya50f10+l'Reduced Travel Time (s)' -BWeSn+t'Local-Regional Shear Phases jsa4.1_example' -K > $psfile
# Reduction velocity
gmt pstext -: -R $projection -F+f12p,Helvetica-Bold,black+a0.+jBL -K -O << END >> $psfile
  -115.000     1.000 Inverse reduction velocity: 31.70 sec/deg
END
gmt psxy -: -R $projection -Wthick,cyan -K -O << END >> $psfile
    -0.100     0.500
    -0.200     1.000
    -0.300     1.500
    -0.400     2.000
    -0.500     2.500
    -0.600     3.000
    -0.700     3.500
    -0.800     4.000
    -0.900     4.500
    -1.000     5.000
    -1.100     5.500
    -1.200     6.000
    -1.300     6.500
    -1.400     7.000
    -1.500     7.500
    -1.600     8.000
    -1.700     8.500
    -1.800     9.000
    -1.900     9.500
    -2.000    10.000
    -2.100    10.500
    -2.200    11.000
    -2.300    11.500
    -2.400    12.000
    -2.500    12.500
    -2.600    13.000
    -2.700    13.500
    -2.800    14.000
    -2.900    14.500
    -3.000    15.000
END
gmt psxy -: -R $projection -Wthick,green -K -O << END >> $psfile
     2.519     1.000
    -0.989     1.500
    -4.498     2.000
    -8.008     2.500
   -11.519     3.000
   -15.032     3.500
   -18.548     4.000
   -22.066     4.500
   -25.587     5.000
   -29.112     5.500
   -32.640     6.000
   -36.173     6.500
   -39.711     7.000
   -43.254     7.500
   -46.803     8.000
   -50.358     8.500
   -53.919     9.000
   -57.487     9.500
   -61.062    10.000
   -64.645    10.500
   -68.236    11.000
   -71.835    11.500
   -75.444    12.000
   -79.062    12.500
   -82.689    13.000
   -86.327    13.500
   -89.975    14.000
   -93.635    14.500
   -97.305    15.000
END
gmt psxy -: -R $projection -Wthick,blue -K -O << END >> $psfile
     0.914     1.000
    -0.540     1.500
    -1.996     2.000
    -3.452     2.500
    -4.909     3.000
    -6.367     3.500
    -7.827     4.000
END
gmt psxy -: -R $projection -Wthick,red -K -O << END >> $psfile
     0.486     0.500
     0.547     1.000
     0.703     1.500
     0.886     2.000
     1.076     2.500
     1.268     3.000
     1.462     3.500
     1.655     4.000
END
gmt psxy -: -R $projection -Sc -Wthin,cyan -K -O << END >> $psfile
END
gmt psxy -: -R $projection -Sx -Wthin,green -K -O << END >> $psfile
   -55.470     7.642     0.200
   -55.365     7.639     0.200
   -57.661     7.869     0.200
   -57.356     7.891     0.200
  -105.408    10.196     0.200
  -108.073    10.248     0.200
  -106.925    10.243     0.200
  -108.019    10.341     0.200
  -109.783    10.428     0.200
   -74.529    10.173     0.200
   -78.206    10.226     0.200
   -76.078    10.222     0.200
   -78.152    10.319     0.200
   -77.886    10.405     0.200
END
gmt psxy -: -R $projection -Sx -Wthin,blue -K -O << END >> $psfile
END
gmt psxy -: -R $projection -Sx -Wthin,red -O << END >> $psfile
END
