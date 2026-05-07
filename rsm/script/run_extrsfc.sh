#!/bin/bash
hxs=0
FEND=120
INCHOUR=1
hr=` expr $hxs + 0 `
while [ $hr -le $FEND ]; do

if [ $hr -lt 10 ] ; then hr=0$hr ; fi
####
#extract sfc variables from pgb files
# usage: extrsfc.sh $hr

if [ $hr = "00" ]; then
 wgrib -d 284 -grib -verf -o extr_sfc.grib r_pgb.f$hr #sfc hgt
 wgrib -d 282 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #mwl u
 wgrib -d 283 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #mwl v
else
 wgrib -d 263 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #sfc pres
 wgrib -d 274 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #sfc CAPE
 wgrib -d 276 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #best LifInd
 wgrib -d 282 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #mwl u
 wgrib -d 283 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #mwl v
 wgrib -d 320 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #APCP
 wgrib -d 325 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #UGRD 10-m
 wgrib -d 326 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #VGRD 10-m
 wgrib -d 327 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #TMP 2-m
 wgrib -d 328 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #SPFH 2-m
 wgrib -d 329 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #TMAX 2-m
 wgrib -d 331 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #TMIN 2-m
 wgrib -d 338 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #HPBL
 wgrib -d 340 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #tcdclls ??? returning albedo
 wgrib -d 299 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #soil T 0-10m
 wgrib -d 300 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #soil T 10-200m
 wgrib -d 297 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #soil W 0-10m
 wgrib -d 298 -grib -verf -append -o extr_sfc.grib r_pgb.f$hr #soil w 10-200m
fi
#####
  hr=`expr $hr + $INCHOUR`
done
grib2ctl.pl -verf extr_sfc.grib > extr_sfc.grib.ctl

exit

