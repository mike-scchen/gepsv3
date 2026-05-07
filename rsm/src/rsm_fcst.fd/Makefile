EXEC = rsm.x

SHELL = /bin/sh
CP = cp
RM = rm -f
CD = cd
CHMOD = chmod

RANLIB = ranlib
AR = ar cq

CPP             =$(C_PREP_FCST)
F90             =$(F_COMP_FCST)
LIBS            =$(F_LIBS_RINP) $(F_LIBS_FCST)
FORT_FLAGS      =$(F_FLAG_FCST) $(F_FLAG)
#W3LB            =$(F_W3LB_FCST)
W3LB            =
INCLUDE_DIR     =-I.

##### rinp_base #####
OBJ_MOD = machine.o  istack.o cnttime.o physcons.o funcphys.o module.o module_progtm.o \
        module_gfsdata.o module_rinpbase.o module_rinpinit.o module_wind_fitch.o \
        module_sw.o module_ras.o module_rascnv.o physics_module.o sigmap_cons.o \
        namelist_soilveg.o \
        radlw_rrtm_param.o radlw_rrtm_datatb.o radlw_rrtm_main.o  radiation_module.o \
        module_mp_wsm3.o phys_mps_wsm5_p.o module_mp_radar.o module_mp_wsm6_fv3.o mod_rndsl.o\
        module_mp_gfdl_fv3_v16.o module_mp_gfdl_fv3_v2.o module_mp_gfdl_fv3_v3.o module_mp_gsfcgce_3ice_nuwrf.o\
        module_bfmicrophysics.o write_to_grads.o \
        mo_constants.o mo_cumulus_flux.o module_SurfaceAlbedo.o sfcsub.o \
        mpi_def.o cmp.comm.o rsm.comm.o
OBJ_FFT99M= fax.o fftrig.o vpassm.o fft99m.o inp_fft99m.o int_fft99m.o fft99a.o fft99b.o

OBJ_RINPBASE= newps.o sg2sg.o trispl.o valts.o setlatlon.o \
	todxdy.o cgtorg.o rtorini.o rclima.o
	
##### rinp_init #####
OBJ_RINPINIT= 
##### rinp_fcst #####

OBJ_RSM = mainrsmf_mrg.o rgetcon_mrg.o gncpus.o rsmini_mpi.o rzerout.o rsmltb.o \
        rloopmp.o rloopsp.o \
        rsmsmf.o rfwdstep.o rsmstep.o rsmdfi.o rdfini.o rwrtout.o rsmsav.o gsm2bgd.o \
        sphpt1.o sphpt2.o sphqk1.o sphqk1mp.o          rpln2i.o ggtorg.o ggtobg.o \
        gtobini.o g2bini.o fplfml.o flpflm.o setbgd.o setg2b.o pint1.o pint3.o getwt.o getwl.o \
        c2bini.o ctobini.o cgtobg.o rsm2bgd.o ll2xy.o xy2ll.o \
        ilagw.o bspline.o tay.o rsetsig.o rsetmap.o rsetsfc.o rfixio.o rsfcmrg.o \
        g2rini_mpi.o c2rini_mpi.o gg2rg_mpi.o cg2rg_mpi.o sread2.o \
        swrite.o swrite2.o getptb.o getful.o gdtocc.o gdtocs.o gdtosc.o cctogd.o cstogd.o \
        sctogd.o cread.o readoz_disprd.o sgl2ful.o ful2sgl.o rwrtsfc.o \
        idsdef.o rbmcm.o ramtm.o rprntb.o ramhm_sp.o rbmcm_sp.o rfidi.o rfidi_sp.o rdeldif.o rsicdif.o rsicdif_sp.o rgetcd.o \
        matinv.o rfilt1.o rfilt2.o rlatbnd.o rltbini.o rupdate.o epslxy.o dellcc.o \
        dellsc.o dellcs.o dz2uv.o uv2dz.o plnini.o sumfcs1.o sumffc1.o sumffs1.o sumggs1.o \
        inpint_plnini.o inpint_rftini.o inpint_gdtocc.o inpint_cctogd.o \
        inpint_ffacosx1.o inpint_ffacosy1.o inpint_ffscosy1.o inpint_ffscosx1.o \
        sumggc1.o sumgcs1.o rmsgtr.o rbartr.o rftini.o \
        ffsncsy1.o ffsncsx1.o ffancsy1.o ffancsx1.o ffscosy1.o ffscosx1.o \
        ffacosy1.o ffacosx1.o \
        ffssiny1.o ffssinx1.o ffasiny1.o ffasinx1.o \
        ffsncsy.o ffsncsx.o ffsncs.o ffancs.o \
        ffacos.o ffasin.o ffscos.o ffssin.o  \
        rimpadj.o rbphys.o extrm2.o setindxoz.o expnd.o\
        ozinterpol.o rsetbgt.o rbgt000.o rbgtlpa.o \
        rbgtsic.o rbgtlpd.o rbgtdif.o rbgtlpb.o rbgtzero.o rbgtini.o rbgtio.o rwrtbgt.o \
        maxmin.o sums2i.o sums2imp.o fl2i.o fl2ic.o fl2imp.o fti_lonf.o gft_lonf.o ggozrm.o \
        sumtop.o sumtopmp.o dzuvle.o dzuvlemp.o dellnp.o pln2i.o gpln2i.o epsilo.o \
        glats.o poly.o cmpind.o transo.o tridi2.o tridin.o tridit.o tridi2t3.o omegast3.o omegast3_sp.o \
        get_prs.o get_phi.o get_prs_sp.o get_phi_sp.o cdefghm.o rloopw.o nfidi.o nsicdif.o hydrow.o nsicini.o newidate.o \
        sreadwav.o swritewav.o writebase.o readbase.o rmpbf2p.o rmpbp2f.o rsm_cpl_init.o rsm_cpl_send2rocn.o rsm_cpl_recv4rocn.o

OBJ_RAD= albaer.o gcljms.o astronomy.o cdate.o compjd.o coszmn.o fcstim.o prtime.o \
        solar.o o3intpnasa.o grrad.o conrad.o crhtab.o rdlwsw.o clo89.o e1e290.o fst88.o \
        lwr88.o e290.o hconst.o spa88.o lwtable.o swr95.o swr95a.o soluv.o solir.o swflux.o \
        flxco2.o aeros.o getaer.o ralsct.o gaeros.o getozn.o ozphys.o \
        getclds.o cldjms.o cldprp.o

OBJ_PHY= progtm.o progt2.o dcyc2t3.o cnvc90.o \
	sfc_diff.o sfc_ocean.o sfc_drv.o sfc_land.o sfc_sice.o sfc_cice.o sfc_diag.o sflx.o set_soilveg.o \
	lrgscl.o moninp.o sascnv.o shalcv.o mstadbt3.o gscond.o precpd.o \
	gscond_nsas.o precpd_nsas.o rgwdps.o mstcnv.o gsmdrive.o \
	micro_init.o sascnv_nsas.o shalcv_nsas.o gscond_sca.o precpd_sca.o \
	mfdeepcnv.o mfshalcnv.o samfdeepcnv.o samfshalcnv.o qsatq.o rayleifr.o \
	moninedmf.o mfpbl.o bl_ysu.o sub_ysu.o lightning_ecn.o \
	cumastr_driv_n.o cuadjtq_n.o cuasc_n.o cutype_n.o cubasmc_n.o cuddraf_n.o cudlfs_n.o \
	cudtdq_n.o cududv_n.o cuentr_n.o cuflx_n.o cuini_n.o cumastr_n.o clima_f2p.o fix_neg_q6.o

OBJ_RAS= setras.o rascnv.o cloud.o crtwrk.o ddrft.o qsatcn.o angrad.o setqrp.o \
        qrabf.o setvtp.o clf.o

OBJ_FUNC= qrpf.o vtpf.o ktsoil.o funcdf.o funckt.o thsat.o twlt.o condense.o deposit.o \
          nTDK_func.o

OBJ_RMPI = dimset.o mpdimset.o mpgp2f.o mpn2nn.o mpnk2nn.o mpnk2nx.o \
mpnl2ny.o mpnn2n.o mpnn2nk.o mpny2nl.o mpsf2p.o mpsp2f.o mpx2nx.o \
spcshfli.o spcshflo.o \
rdimset.o rmpbcasti.o rmpnk2nn.o rmpsynrow.o rmpbcastc.o \
equdis.o equdiv.o rmpbcastr.o rmpbcastl.o rmpgf2p.o rmpnk2nx.o rmpny2nl.o  \
rmpgp2f.o rmpnl2ny.o \
rmpfine.o rmpinit.o \
rmpabort.o rmpnn2nk.o rmpsynall.o \
rmpnx2nk.o rmpsyncol.o \
rmpsf2p.o rmpsp2f.o rmpgetspd.o \
rmpdimset.o rmpnx2x.o rmpn2nn.o rmpnn2n.o rmpx2nx.o rmpkencollect.o \
recv_data.o recv_idate.o RSM_RECVA2Rdata.o read_data.o read_idate.o mpi_more.o

OBJ_RKEN = rkenpre.o rgetij.o rkenput.o rkenputr.o rkenwrt.o rsmkenprof.o

OBJ_W3 = w3movdat.o w3fs26.o iw3jdn.o w3reddat.o
OBJ_NDSL= rndslfv_monoadvh2.o rndslfv_monoadvv.o rndslfv_monoadvv_nonhyd.o \
					rndslfv_monoadvv_tau.o rndslfv_pack.o rndslfv_para.o dediffusion.o \
					rloopmp_euNDSL.o rfidi_eundsl.o rfidi_sp_eundsl.o nfidi_eundsl.o rdeldif_twice.o

OBJS = ${OBJ_MOD} ${OBJ_FFT99M} ${OBJ_RINPBASE} ${OBJ_RINPINIT} ${OBJ_RSM} ${OBJ_RAD} \
	     ${OBJ_PHY} ${OBJ_RAS} ${OBJ_FUNC} \
       ${OBJ_RMPI} ${OBJ_RKEN} ${OBJ_W3} ${OBJ_NDSL}

SRCS = ${OBJS:.o=.F}
#
# *****************************************************************
#
all: model-rmpi

lib_4cpl_v1: librsm4cpl_v1.a
lib_4cpl: librsm4cpl.a
lib_2cpl: librsm2cpl.a
lib: librsm.a

librsm4cpl_v1.a: ${OBJS}
	@echo
	@echo '==== LINK ===================================================='
	$(RM) $(@)
	$(AR) $(@) $(OBJS)
	$(RANLIB) $(@)
librsm4cpl.a: ${OBJS}
	@echo
	@echo '==== LINK ===================================================='
	$(RM) $(@)
	$(AR) $(@) $(OBJS)
	$(RANLIB) $(@)
librsm2cpl.a: ${OBJS}
	@echo
	@echo '==== LINK ===================================================='
	$(RM) $(@)
	$(AR) $(@) $(OBJS)
	$(RANLIB) $(@)
librsm.a: ${OBJS}
	@echo
	@echo '==== LINK ===================================================='
	$(RM) $(@)
	$(AR) $(@) $(OBJS)
	$(RANLIB) $(@)


model-rmpi: ${OBJS}
	${F90} ${FORT_FLAGS} -o ${EXEC} ${OBJS} ${W3LB} ${LIBS}

clean:
	rm -f *.o *.x *.$(F_SUBF) *.mod core *.lst 

.SUFFIXES : .o .$(F_SUBF) .F .f

.F.$(F_SUBF) : *.h
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed -e '/^\/\*/,/\*\/$$/d' -e '/^ *$$/d' $*.i >$*.$(F_SUBF) ; rm $*.i 

.F.o : *.h
	${CPP} -P ${INCLUDE_DIR} $*.F >$*.i
	sed -e '/^\/\*/,/\*\/$$/d' -e '/^ *$$/d' -e 's/MPI_COMM_WORLD/mpi_comm_rsm/Ig' -e 's/use machine/use rmachine/Ig' -e 's/use physcons/use rphyscons/Ig' -e 's/use funcphys/use rfuncphys/Ig' $*.i >$*.$(F_SUBF) ; rm $*.i
	${F90} ${FORT_FLAGS} -c $*.$(F_SUBF)

.$(F_SUBF).o : *.h
	${F90} ${FORT_FLAGS} -c $*.$(F_SUBF)

.f.o : *.h
	${F90} ${FORT_FLAGS} -c $*.f

#
# *****************************************************************

