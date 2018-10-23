#!/bin/bash
gmt gmtset PS_MEDIA letter
gmt gmtset PROJ_LENGTH_UNIT cm
gmt gmtset FORMAT_GEO_MAP D
gmt gmtset PS_PAGE_ORIENTATION portrait
psfile=jsa4.1_example_base.ps
projection=-JM16.0c+
region=-R-148.4/-147.9/-7.9/-7.2
gmt psbasemap $region $projection -Bxa0.2f0.1 -Bya0.2f0.1 -BWeSn+t'Base Map jsa4.1_example' -K > $psfile
gmt pscoast $projection $region -Df -Ia,blue -Wblue -O -K >> $psfile
# Event numbers
gmt psxy -: -R $projection -Sl8p+t"1" -G0 -O -K << END >> $psfile
    -7.385  -148.341
END
gmt psxy -: -R $projection -Sl8p+t"2" -G0 -O -K << END >> $psfile
    -7.367  -148.298
END
gmt psxy -: -R $projection -Sl8p+t"3" -G0 -O -K << END >> $psfile
    -7.328  -148.349
END
gmt psxy -: -R $projection -Sl8p+t"4" -G0 -O -K << END >> $psfile
    -7.354  -148.321
END
# Change in location, calibration shift if available
gmt psxy -: -R $projection -Wthin,green -O -K << END >> $psfile
    -7.500  -148.262
    -7.385  -148.341
END
gmt psxy -: -R $projection -Wthin,black -O -K << END >> $psfile
    -7.500  -148.262
    -7.385  -148.341
END
gmt psxy -: -R $projection -Wthin,green -O -K << END >> $psfile
    -7.511  -148.225
    -7.367  -148.298
END
gmt psxy -: -R $projection -Wthin,black -O -K << END >> $psfile
    -7.511  -148.225
    -7.367  -148.298
END
gmt psxy -: -R $projection -Wthin,green -O -K << END >> $psfile
    -7.625  -148.060
    -7.328  -148.349
END
gmt psxy -: -R $projection -Wthin,black -O -K << END >> $psfile
    -7.625  -148.060
    -7.328  -148.349
END
gmt psxy -: -R $projection -Wthin,green -O -K << END >> $psfile
    -7.750  -148.040
    -7.354  -148.321
END
gmt psxy -: -R $projection -Wthin,black -O -K << END >> $psfile
    -7.750  -148.040
    -7.354  -148.321
END
# Confidence ellipses for relative location
gmt psxy -: -R $projection -SE -Wthick -O -K << END >> $psfile
    -7.385  -148.341    12.616     4.509    10.388
    -7.367  -148.298    18.967     5.297    12.640
    -7.328  -148.349    12.645     4.628    10.588
    -7.354  -148.321    12.304     4.658    12.606
END
# Circle of radius 5 km for reference
gmt psxy -: -R $projection -SE -Wthick,red -O -K << END >> $psfile
    -7.796  -148.395     0.000    10.000    10.000
END
gmt psxy -: -R $projection -Sx -Wthinnest,blue -O -K << END >> $psfile
    -7.796  -148.395     0.070
END
gmt pstext -: -R $projection -F+f6p,Helvetica-Bold,red+a0.+jBC -Gwhite -O << END >> $psfile
    -7.841  -148.395   5 km
END
