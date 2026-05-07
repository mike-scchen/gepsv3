#ifdef MP
#define LNWAVS lnwavp
#else
#define LNWAVS lnwav
#endif
      common/rcomio/  delx, dely, runid, usrid,                          &    
     &                snnp1(LNWAVS),rnnp1(LNWAVS),rnnp1max,		 &
     &                epsx(LNWAVS),epsy(LNWAVS),epsxmax,epsymax,	 &
     &                dpsx(LNWAVS),dpsy(LNWAVS),			 &
     &                deltim,dt2,dtltb,zhour,phour,fhour,odt_cpl,        &
     &                rcl,sl1,dk,tk
      common/rcomioi/ iimprlx,isemimp,idmpjet,imdlphy,			 &
     &                ifin,icen,igen,icen2,ienst,iensi,ndigyr,		 &
     &                kdt,limlow,numsum,nummax,ncldb1,irclim,kdt_cpl
      common/rcomioi/ nmtnv,nctune,nco2,nozon,				 &
     &                nrslm,nrorg,                     &
#ifdef  A
     &                nrbgt,nrbgt1,nrbgt2,				 &
#endif
     &                nrinit,nrflxi,nrkenp,				 &
     &                nrsmi0,nrsmi1,nrsmi2,nrflip,nb, &
     &                nrsmop,nrsfcp,nrflxp,				 &
     &                nrsmo1,nrsmo2,nrflop,        &
     &                nrsms1,nrsms2,nrflop8
      common/rcomio/  rlxmsec,rlxhsec,filta,filtb,seminpa,		 &
     &                difmsec,difhsec,rsfcsec,restrhr
      common/rcomio/  fhseg,fhout,fhswr,fhlwr,fhbas,fhdfi,fhken
      common/rcomioi/ nsseg,nsout,nsswr,nslwr,nsbas,nsdfi,nsclm
      logical         lfnhr,lsfwd,lsftr,lssav,lsout,lsimp,lsphy,         &  
     &                lscca,lsswr,lslwr,lsbgt,lsken,lsdmp,lsfcmrg,       &
     &                lrestart,lrclm,lsclm
      common/rcomiol/ lfnhr,lsfwd,lsftr,lssav,lsout,lsimp,lsphy,         &  
     &                lscca,lsswr,lslwr,lsbgt,lsken,lsdmp,lsfcmrg,       &
     &                lrestart,lrclm,lsclm
