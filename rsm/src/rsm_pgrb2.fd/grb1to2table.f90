subroutine grb1to2table(IPU,p0,p1,p2)
integer::ip0(255),ip1(255),ip2(255)
integer::p0,p1,p2
!===========
data ip0(  1),ip1(  1),ip2(  1) /0,3, 0/  !PRES
data ip0(  2),ip1(  2),ip2(  2) /0,3, 1/  !PRMSL
data ip0(  3),ip1(  3),ip2(  3) /0,3, 2/  !PTEND 
data ip0(  7),ip1(  7),ip2(  7) /0,3, 5/  !HGT
data ip0( 11),ip1( 11),ip2( 11) /0,0, 0/  !TMP
data ip0( 15),ip1( 15),ip2( 15) /0,0, 4/  !TMAX
data ip0( 16),ip1( 16),ip2( 16) /0,0, 5/  !TMIN
data ip0( 17),ip1( 17),ip2( 17) /0,0, 6/  !DPT
data ip0( 33),ip1( 33),ip2( 33) /0,2, 2/  !UGRD
data ip0( 34),ip1( 34),ip2( 34) /0,2, 3/  !VGRD
data ip0( 39),ip1( 39),ip2( 39) /0,2, 8/  !VVEL
data ip0( 40),ip1( 40),ip2( 40) /0,2, 9/  !DZDT
data ip0( 41),ip1( 41),ip2( 41) /0,2,10/  !ABSV
data ip0( 51),ip1( 51),ip2( 51) /0,1, 0/  !SPFH
data ip0( 52),ip1( 52),ip2( 52) /0,1, 1/  !RH
data ip0( 54),ip1( 54),ip2( 54) /0,1, 3/  !PWAT
data ip0( 61),ip1( 61),ip2( 61) /0,1, 8/  !APCP
data ip0( 63),ip1( 63),ip2( 63) /0,1,10/  !ACPCP
data ip0( 65),ip1( 65),ip2( 65) /0,1,13/  !WEASD
data ip0( 66),ip1( 66),ip2( 66) /0,1,11/  !SNOD
data ip0( 71),ip1( 71),ip2( 71) /0,6, 1/  !TCDC
data ip0( 81),ip1( 81),ip2( 81) /2,0, 0/  !LAND
data ip0( 84),ip1( 84),ip2( 84) /0,19,1/  !ALBEDO
data ip0( 90),ip1( 90),ip2( 90) /2 ,0,5/  !WATR water runoff
data ip0( 91),ip1( 91),ip2( 91) /10,2,0/  !ICEC
data ip0( 85),ip1( 85),ip2( 85) /2,0, 3/  !TSOIL
data ip0( 86),ip1( 86),ip2( 86) /2,3,19/  !SOILMOI
data ip0(121),ip1(121),ip2(121) /0,0,10/  !LHTFL
data ip0(122),ip1(122),ip2(122) /0,0,11/  !SHTFL
data ip0(124),ip1(124),ip2(124) /0,2,17/  !UFLX
data ip0(125),ip1(125),ip2(125) /0,2,18/  !VFLX
data ip0(131),ip1(131),ip2(131) /0,7,10/  !LFTX
data ip0(132),ip1(132),ip2(132) /0,7,11/  !4LFTX
data ip0(136),ip1(136),ip2(136) /0,2,25/  !VWSH
data ip0(140),ip1(140),ip2(140) /0,1,192/ !CRAIN
data ip0(141),ip1(141),ip2(141) /0,1,193/ !CFRZR
data ip0(142),ip1(142),ip2(142) /0,1,194/ !CICEP
data ip0(143),ip1(143),ip2(143) /0,1,195/ !CSNOW
data ip0(144),ip1(144),ip2(144) /2,0, 9/  !SOILW
data ip0(145),ip1(145),ip2(145) /0,1,41/  !PEVPR Potential Evaporation Rate
data ip0(146),ip1(146),ip2(146) /0,6,15/  !CWORK
data ip0(147),ip1(147),ip2(147) /0,2,210/ !U-GWD
data ip0(148),ip1(148),ip2(148) /0,2,211/ !V-GWD
data ip0(153),ip1(153),ip2(153) /0,1,22/  !CLWMR
data ip0(154),ip1(154),ip2(154) /0,14,1/  !O3MR
data ip0(155),ip1(155),ip2(155) /2,0,10/  !GFLUX
data ip0(156),ip1(156),ip2(156) /0,7, 7/  !CIN
data ip0(157),ip1(157),ip2(157) /0,7, 6/  !CAPE
data ip0(170),ip1(170),ip2(170) /0,1,24/  !RWMR
data ip0(171),ip1(171),ip2(171) /0,1,25/  !SNMR
data ip0(178),ip1(178),ip2(178) /0,1,23/  !ICMR
data ip0(179),ip1(179),ip2(179) /0,1,32/  !GRLE
data ip0(180),ip1(180),ip2(180) /0,2,22/  !GUST
data ip0(204),ip1(204),ip2(204) /0,4, 7/  !DSWRF
data ip0(205),ip1(205),ip2(205) /0,5, 3/  !DLWRF
data ip0(211),ip1(211),ip2(211) /0,4, 8/  !USWRF
data ip0(212),ip1(212),ip2(212) /0,5, 4/  !ULWRF
data ip0(221),ip1(221),ip2(221) /0,3,18/  !HPBL 
data ip0(222),ip1(222),ip2(222) /0,3,15/  !5WAVH
data ip0(235),ip1(235),ip2(235) /1,0, 6/  !SSRUN
data ip0(242),ip1(242),ip2(242) /2,0,224/ !WVUFLX
data ip0(243),ip1(243),ip2(243) /2,0,225/ !WVVFLX

p0=ip0(IPU)
p1=ip1(IPU)
p2=ip2(IPU)
!print*,p0,p1,p2
endsubroutine
