#!/bin/bash
kshpath=`pwd `
. ${kshpath}/00_setDate.sh


if [ $1 = 639  ] ; then
  MODEL_BASIC='nco=640,'
  MODLST_RES='dt=225., hfilt=1., cgw=4.2e-5,'
  DMSFLAG=GJ
  NPEX=4
  NPEY=384
elif [ $1 = 383  ] ; then
  MODEL_BASIC='nco=384,'
  MODLST_RES="dt=600., hfilt=1., cgwd=2.4, cmbk=0.6, tofd=t, ksgeo=1, outdms=0, outgrb2=1,
              domfc=1080., out_green=t, otgreen=3., nmmiph=12, nmgwcv=1, two_loop=t,
              nmcup=7, nmshl=4,factop=80.,mass_dp=t,itter=1,alpha=0.65,"
  DMSFLAG=GI
  NPEX=4
  NPEY=64
fi

#-- write out file list
cat > ${GFSWRK}/filist << EOF
&filst
 ifilin=ANADMS,
 ifilout=FCSTDMS,
 cwbout='${GFSWRK}/tmpdir/gfs_cwbout',
 bckfile=BCKOPS,
 phyout='${GFSWRK}/tmpdir/gfs_phyout',
 namlsts='${GFSWRK}/namlsts',
 crdate='${GFSWRK}/crdate',
 ocards='${GFSWRK}/ocards',
 cntrl='${GFSWRK}/gfsctl',
 ifilout_grb='${GFSOUT}',
&end
EOF

#-- write out running date tag
echo $dtg > ${GFSWRK}/crdate

#-- write out gfsctl
cat ${GFSDIR2}/etc_4cpl/gfsctl > ${GFSWRK}/gfsctl

#-- write out ocards
cat ${GFSDIR2}/etc_4cpl/ocards > ${GFSWRK}/ocards

#-- write out namlsts
cat > ${GFSWRK}/namlsts << EOF
&model_param
 ${MODEL_BASIC}
 lev=72,
 ncld=7,
 octahedral=t,
 nout=90000,
 io_quilting=f,
 npex=${NPEX},
 npey=${NPEY},
&end

&modlst
if [ ${struc} = 4cpl  ] ; then
 outrsm=true, rsmoutinv=6, rlon1=100., rlon2=150., rlat1=5., rlat2=40., rgrdsz=0.25,
else
 outrsm=false, rsmoutinv=6, rlon1=100., rlon2=150., rlat1=5., rlat2=40., rgrdsz=0.25,
fi
 taui=0.0, taue=120.0, tauo=1.0, taup=6.0, taureg=0.,
 cstar=f, update=t, lsimpl=t,
 ksgeo=1, yesdia=t,
 dopbl=t, docup=t, dorad=t, dolsp=t, doshl=t, dodry=f,
 dograv=t, docgrav=t,
 donnmi=t,
 dosppt=t, doskeb=t, dospptout=f,
 doshum=f,
 cutfreq=3, nnmivm=3,
 doincr=f,
 hdiff=t, frad=1.0, ldiag=0,
 idg=40, jdg=108,
 itypbl=0, numreduce=8, ptmeans=800.,
 irad=2, nmland=2,
 nmcup=6, nmshl=3, nmpbl=4, nmmiph=2,
 nmgwor=2, nmgwcv=1,
 ktcup=20,
 mtnvar=14, doo3l=t,
 ioutsigr=1,
 ggdef='${DMSFLAG}0G', gmdef='${DMSFLAG}MG',
 domfc=1080., out_green=t, otgreen=3., out_hp=f,
 ndsladvh2=f,
 isot=1, ivegsrc=1, cgwd=2.40, cmbk=0.60,
 spl1=5., spl2=100., tofd=t
 ${MODLST_RES}
 &end

&typ
 write_mem=0,
 trk_intv=6,
&end

&stochy_physics
 ncep_seeds = true,
 use_zmtnblck = true,
 sppt_logit = true,
 sppt_sigtop1 = 0.1,
 sppt_sigtop2 = 0.025,
 sppt_sfclimit = true,
 sppt_sigbot1 = 0.975,
 sppt_sigbot2 = 0.9,
 sppt = 0.70,0.4,0.10,0.08,0.04
 sppt_seed = -999,-999,-999,-999,-999
 sppt_decort = 2.16E4,2.592E5,2.592E6,7.776E6,3.1536E7
 sppt_lscale = 500.E3,1000.E3,2000.E3,2000.E3,2000.E3
 shum = 0.18,-999,-999,-999,-999
 shum_seed = -999,-999,-999,-999,-999
 shum_decort = 2.16E4,1.728E5,2.592E6,7.776E6,3.1536E7
 shum_lscale = 500.E3,1000.E3,2000.E3,2000.E3,2000.E3
 shum_sigefold = 0.2,
 skeb_sigtop1 = 0.1,
 skeb_sigtop2 = 0.025,
 skeb_sigbot1 = 0.975,
 skeb_sigbot2 = 0.9,
 skeb_vdof = 5,
 skebnorm = 1,
 skebfilt = 12,
 skeb = 1.E12,-999,-999,-999,-999
 skeb_seed = -999,-999,-999,-999,-999
 skeb_decort = 2.16E4,2.592E5,2.592E6,7.776E6,3.1536E7
 skeb_lscale = 500.E3,1000.E3,2000.E3,2000.E3,2000.E3
 ssst = 0.80,0.4,0.10,0.08,0.04
 ssst_seed = -999,-999,-999,-999,-999
 ssst_decort = 2.16E4,2.592E5,2.592E6,7.776E6,3.1536E7
 ssst_lscale = 500.E3,1000.E3,2000.E3,2000.E3,2000.E3
&end

&grb_conf
 grbmem=002,
 grbnumm=${member},
&end

&gce_3ice
 SL_sedi=false, sat_predict=true, new_saturation=true,
 use_cpm=false, use_declination=false
&end

EOF

ln -fs ${GFSDIR}/fix/* ${GFSWRK}/
#ln -fs ${GFSWRK}/global_o3prdlos.f77 ${GFSWRK}/fort.28
ln -fs ${GFSWRK}/ozprdlos_2015_new_sbuvO3_tclm15_nuchem.f77  ${GFSWRK}/global_o3prdlos
ln -fs ${GFSWRK}/global_o3clim.txt   ${GFSWRK}/fort.48
ln -fs ${GFSWRK}/global_climaeropac_global.txt ${GFSWRK}/aerosol.dat
ln -fs ${GFSWRK}/global_sfc_emissivity_idx.txt ${GFSWRK}/sfc_emissivity_idx.txt

