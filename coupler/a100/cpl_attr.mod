V34 :0x34 cpl_attr
12 cpl_attr.F90 S624 0
10/27/2025  15:42:24
use cpl_rank public 0 direct
use m_zeit private
use m_list private
use m_router private
use m_attrvect private
enduse
D 196 26 1033 1104 1030 7
D 329 23 6 1 10 382 0 0 0 0 0
 10 381 11 10 381 382
D 332 23 6 1 11 382 0 0 0 0 0
 0 382 11 11 382 382
D 337 20 168
D 339 23 337 1 10 382 0 0 0 0 0
 10 381 11 10 381 382
D 342 23 337 1 11 382 0 0 0 0 0
 0 382 11 11 382 382
D 502 20 168
D 540 26 1858 2208 1857 7
D 660 20 726
S 624 24 0 0 0 10 1 0 5012 10005 8000 A 0 0 0 0 B 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 cpl_attr
S 626 23 0 0 0 10 1030 624 5032 4 0 A 0 0 0 0 B 400000 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 attrvect
S 627 19 0 0 0 10 1 624 5041 4 0 A 0 0 0 0 B 400000 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 192 3 0 0 0 0 0 624 0 0 0 0 attrvect_init
O 627 3 1161 1155 1150
S 629 19 0 0 0 10 1 624 5060 4 0 A 0 0 0 0 B 400000 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 188 2 0 0 0 0 0 624 0 0 0 0 attrvect_importrattr
O 629 2 1291 1284
S 631 19 0 0 0 6 1 624 5093 4 0 A 0 0 0 0 B 400000 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 185 1 0 0 0 0 0 624 0 0 0 0 mct_avect_indexia
O 631 1 1208
S 633 19 0 0 0 6 1 624 5119 4 0 A 0 0 0 0 B 400000 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 183 1 0 0 0 0 0 624 0 0 0 0 mct_avect_indexra
O 633 1 1215
S 636 23 0 0 0 10 1857 624 5154 4 0 A 0 0 0 0 B 400000 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 router
S 637 19 0 0 0 10 1 624 5161 4 0 A 0 0 0 0 B 400000 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 270 2 0 0 0 0 0 624 0 0 0 0 router_init
O 637 2 1971 1964
R 841 26 24 m_list =
S 1009 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 1010 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 1030 25 5 m_attrvect attrvect
R 1033 5 8 m_attrvect ilist attrvect
R 1034 5 9 m_attrvect rlist attrvect
R 1037 5 12 m_attrvect iattr attrvect
R 1038 5 13 m_attrvect iattr$sd attrvect
R 1039 5 14 m_attrvect iattr$p attrvect
R 1040 5 15 m_attrvect iattr$o attrvect
R 1044 5 19 m_attrvect rattr attrvect
R 1045 5 20 m_attrvect rattr$sd attrvect
R 1046 5 21 m_attrvect rattr$p attrvect
R 1047 5 22 m_attrvect rattr$o attrvect
R 1150 14 125 m_attrvect init_
R 1155 14 130 m_attrvect initv_
R 1161 14 136 m_attrvect initl_
R 1208 14 183 m_attrvect indexia_
R 1215 14 190 m_attrvect indexra_
R 1284 14 259 m_attrvect importrattrsp_
R 1291 14 266 m_attrvect importrattrdp_
S 1404 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 1405 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 16 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 1406 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 1409 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 1410 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 1413 3 0 0 0 502 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 8220 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 8 5b 4d 57 54 49 4d 45 5d
S 1414 3 0 0 0 502 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 8229 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 8 5b 58 57 54 49 4d 45 5d
S 1415 3 0 0 0 502 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 8238 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 8 5b 50 55 54 49 4d 45 5d
S 1416 3 0 0 0 502 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 8247 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 8 5b 50 53 54 49 4d 45 5d
S 1417 3 0 0 0 502 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 8256 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 8 5b 43 55 54 49 4d 45 5d
S 1418 3 0 0 0 502 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 8265 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 8 5b 43 53 54 49 4d 45 5d
R 1448 7 21 m_zeit masks$ac
R 1451 7 24 m_zeit header$ac
R 1857 25 18 m_router router
R 1858 5 19 m_router comp1id router
R 1859 5 20 m_router comp2id router
R 1860 5 21 m_router nprocs router
R 1861 5 22 m_router maxsize router
R 1862 5 23 m_router lavsize router
R 1863 5 24 m_router numiatt router
R 1864 5 25 m_router numratt router
R 1866 5 27 m_router pe_list router
R 1867 5 28 m_router pe_list$sd router
R 1868 5 29 m_router pe_list$p router
R 1869 5 30 m_router pe_list$o router
R 1872 5 33 m_router pe_list_loc router
R 1873 5 34 m_router pe_list_loc$sd router
R 1874 5 35 m_router pe_list_loc$p router
R 1875 5 36 m_router pe_list_loc$o router
R 1877 5 38 m_router mpicomm router
R 1879 5 40 m_router num_segs router
R 1880 5 41 m_router num_segs$sd router
R 1881 5 42 m_router num_segs$p router
R 1882 5 43 m_router num_segs$o router
R 1885 5 46 m_router locsize router
R 1886 5 47 m_router locsize$sd router
R 1887 5 48 m_router locsize$p router
R 1888 5 49 m_router locsize$o router
R 1891 5 52 m_router permarr router
R 1892 5 53 m_router permarr$sd router
R 1893 5 54 m_router permarr$p router
R 1894 5 55 m_router permarr$o router
R 1898 5 59 m_router seg_starts router
R 1899 5 60 m_router seg_starts$sd router
R 1900 5 61 m_router seg_starts$p router
R 1901 5 62 m_router seg_starts$o router
R 1905 5 66 m_router seg_lengths router
R 1906 5 67 m_router seg_lengths$sd router
R 1907 5 68 m_router seg_lengths$p router
R 1908 5 69 m_router seg_lengths$o router
R 1911 5 72 m_router rp1 router
R 1912 5 73 m_router rp1$sd router
R 1913 5 74 m_router rp1$p router
R 1914 5 75 m_router rp1$o router
R 1917 5 78 m_router ip1 router
R 1918 5 79 m_router ip1$sd router
R 1919 5 80 m_router ip1$p router
R 1920 5 81 m_router ip1$o router
R 1923 5 84 m_router ireqs router
R 1924 5 85 m_router ireqs$sd router
R 1925 5 86 m_router ireqs$p router
R 1926 5 87 m_router ireqs$o router
R 1928 5 89 m_router rreqs router
R 1930 5 91 m_router rreqs$sd router
R 1931 5 92 m_router rreqs$p router
R 1932 5 93 m_router rreqs$o router
R 1936 5 97 m_router istatus router
R 1937 5 98 m_router istatus$sd router
R 1938 5 99 m_router istatus$p router
R 1939 5 100 m_router istatus$o router
R 1941 5 102 m_router rstatus router
R 1944 5 105 m_router rstatus$sd router
R 1945 5 106 m_router rstatus$p router
R 1946 5 107 m_router rstatus$o router
R 1964 14 125 m_router initd_
R 1971 14 132 m_router initp_
S 1996 6 4 0 0 196 1997 624 10787 4 8 A 0 0 0 0 B 0 12 0 0 0 0 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfs_export_attr
S 1997 6 4 0 0 196 1998 624 10803 4 8 A 0 0 0 0 B 0 13 0 0 0 1104 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfs_send_gocn_attr
S 1998 6 4 0 0 196 1999 624 10822 4 8 A 0 0 0 0 B 0 14 0 0 0 2208 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfs_recv_gocn_attr
S 1999 6 4 0 0 196 2000 624 10841 4 8 A 0 0 0 0 B 0 15 0 0 0 3312 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfs_nested_attr
S 2000 6 4 0 0 196 2001 624 10857 4 8 A 0 0 0 0 B 0 16 0 0 0 4416 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfs_send_rsm_attr
S 2001 6 4 0 0 196 2002 624 10875 4 8 A 0 0 0 0 B 0 18 0 0 0 5520 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rsm_export_attr
S 2002 6 4 0 0 196 2003 624 10891 4 8 A 0 0 0 0 B 0 19 0 0 0 6624 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rsm_send_rocn_attr
S 2003 6 4 0 0 196 2004 624 10910 4 8 A 0 0 0 0 B 0 20 0 0 0 7728 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rsm_recv_rocn_attr
S 2004 6 4 0 0 196 2005 624 10929 4 8 A 0 0 0 0 B 0 21 0 0 0 8832 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rsm_nested_attr
S 2005 6 4 0 0 196 2006 624 10945 4 8 A 0 0 0 0 B 0 22 0 0 0 9936 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rsm_recv_gfs_attr
S 2006 6 4 0 0 196 2007 624 10963 4 8 A 0 0 0 0 B 0 24 0 0 0 11040 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocn_export_attr
S 2007 6 4 0 0 196 2008 624 10980 4 8 A 0 0 0 0 B 0 25 0 0 0 12144 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocn_send_gfs_attr
S 2008 6 4 0 0 196 2009 624 10999 4 8 A 0 0 0 0 B 0 26 0 0 0 13248 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocn_recv_gfs_attr
S 2009 6 4 0 0 196 2010 624 11018 4 8 A 0 0 0 0 B 0 27 0 0 0 14352 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocn_nested_attr
S 2010 6 4 0 0 196 2011 624 11035 4 8 A 0 0 0 0 B 0 28 0 0 0 15456 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocn_send_rocn_attr
S 2011 6 4 0 0 196 2012 624 11055 4 8 A 0 0 0 0 B 0 30 0 0 0 16560 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rocn_export_attr
S 2012 6 4 0 0 196 2013 624 11072 4 8 A 0 0 0 0 B 0 31 0 0 0 17664 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rocn_send_rsm_attr
S 2013 6 4 0 0 196 2014 624 11091 4 8 A 0 0 0 0 B 0 32 0 0 0 18768 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rocn_recv_rsm_attr
S 2014 6 4 0 0 196 2015 624 11110 4 8 A 0 0 0 0 B 0 33 0 0 0 19872 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rocn_nested_attr
S 2015 6 4 0 0 196 2016 624 11127 4 8 A 0 0 0 0 B 0 34 0 0 0 20976 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rocn_recv_gocn_attr
S 2016 6 4 0 0 540 2017 624 11147 4 8 A 0 0 0 0 B 0 36 0 0 0 22080 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfs_send_gocn_rout
S 2017 6 4 0 0 540 2018 624 11166 4 8 A 0 0 0 0 B 0 37 0 0 0 24288 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfs_recv_gocn_rout
S 2018 6 4 0 0 540 2019 624 11185 4 8 A 0 0 0 0 B 0 38 0 0 0 26496 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfs_send_rsm_rout
S 2019 6 4 0 0 540 2020 624 11203 4 8 A 0 0 0 0 B 0 40 0 0 0 28704 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rsm_send_rocn_rout
S 2020 6 4 0 0 540 2021 624 11222 4 8 A 0 0 0 0 B 0 41 0 0 0 30912 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rsm_recv_rocn_rout
S 2021 6 4 0 0 540 2022 624 11241 4 8 A 0 0 0 0 B 0 42 0 0 0 33120 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rsm_recv_gfs_rout
S 2022 6 4 0 0 540 2023 624 11259 4 8 A 0 0 0 0 B 0 44 0 0 0 35328 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocn_send_gfs_rout
S 2023 6 4 0 0 540 2024 624 11278 4 8 A 0 0 0 0 B 0 45 0 0 0 37536 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocn_recv_gfs_rout
S 2024 6 4 0 0 540 2025 624 11297 4 8 A 0 0 0 0 B 0 46 0 0 0 39744 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocn_send_rocn_rout
S 2025 6 4 0 0 540 2026 624 11317 4 8 A 0 0 0 0 B 0 48 0 0 0 41952 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rocn_send_rsm_rout
S 2026 6 4 0 0 540 2027 624 11336 4 8 A 0 0 0 0 B 0 49 0 0 0 44160 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rocn_recv_rsm_rout
S 2027 6 4 0 0 540 1 624 11355 4 8 A 0 0 0 0 B 0 50 0 0 0 46368 0 0 0 0 0 0 2061 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rocn_recv_gocn_rout
S 2030 3 0 0 0 6 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 256 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 2031 6 4 0 0 660 2034 624 11428 80000c 8 A 0 0 0 0 B 0 52 0 0 0 0 0 0 0 0 0 0 2062 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfsexport_rlist
S 2034 6 4 0 0 660 2037 624 11494 80000c 8 A 0 0 0 0 B 0 53 0 0 0 256 0 0 0 0 0 0 2062 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rsmexport_rlist
S 2037 6 4 0 0 660 2040 624 11537 80000c 8 A 0 0 0 0 B 0 54 0 0 0 512 0 0 0 0 0 0 2062 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocnexport_rlist
S 2040 6 4 0 0 660 2042 624 11566 80000c 8 A 0 0 0 0 B 0 55 0 0 0 768 0 0 0 0 0 0 2062 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rocnexport_rlist
S 2042 6 4 0 0 660 2043 624 11597 80000c 8 A 0 0 0 0 B 0 57 0 0 0 1024 0 0 0 0 0 0 2062 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfsnested_rlist
S 2043 6 4 0 0 660 2046 624 11613 80000c 8 A 0 0 0 0 B 0 58 0 0 0 1280 0 0 0 0 0 0 2062 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rsmnested_rlist
S 2046 6 4 0 0 660 2047 624 11653 80000c 8 A 0 0 0 0 B 0 59 0 0 0 1536 0 0 0 0 0 0 2062 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocnnested_rlist
S 2047 6 4 0 0 660 2048 624 11670 80000c 8 A 0 0 0 0 B 0 60 0 0 0 1792 0 0 0 0 0 0 2062 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rocnnested_rlist
S 2048 6 4 0 0 660 2049 624 11687 80000c 8 A 0 0 0 0 B 0 63 0 0 0 2048 0 0 0 0 0 0 2062 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfs2gocn_rlist
S 2049 6 4 0 0 660 2050 624 11702 80000c 8 A 0 0 0 0 B 0 64 0 0 0 2304 0 0 0 0 0 0 2062 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocn2gfs_rlist
S 2050 6 4 0 0 660 2051 624 11717 80000c 8 A 0 0 0 0 B 0 65 0 0 0 2560 0 0 0 0 0 0 2062 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfs2rsm_rlist
S 2051 6 4 0 0 660 2052 624 11731 80000c 8 A 0 0 0 0 B 0 66 0 0 0 2816 0 0 0 0 0 0 2062 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocn2rocn_rlist
S 2052 6 4 0 0 660 2053 624 11747 80000c 8 A 0 0 0 0 B 0 67 0 0 0 3072 0 0 0 0 0 0 2062 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rsm2rocn_rlist
S 2053 6 4 0 0 660 1 624 11762 80000c 8 A 0 0 0 0 B 0 68 0 0 0 3328 0 0 0 0 0 0 2062 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rocn2rsm_rlist
S 2060 19 0 0 0 6 1 624 5055 4 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 277 5 0 0 0 0 0 624 0 0 0 0 init
O 2060 5 1971 1964 1161 1155 1150
S 2061 11 0 0 0 10 1995 624 13319 40800000 805000 A 0 0 0 0 B 0 69 0 0 0 48576 0 0 1996 2027 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _cpl_attr$0
S 2062 11 0 0 0 10 2061 624 13331 40800000 805000 A 0 0 0 0 B 0 69 0 0 0 3584 0 0 2031 2053 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _cpl_attr$9
S 2063 23 5 0 0 0 2067 624 13343 0 0 A 0 0 0 0 B 0 71 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 cpl_attr_init
S 2064 1 3 1 0 6 1 2063 13357 4 3000 A 0 0 0 0 B 0 71 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 id_src
S 2065 1 3 1 0 6 1 2063 13364 4 3000 A 0 0 0 0 B 0 71 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 id_dst
S 2066 1 3 1 0 6 1 2063 13371 4 3000 A 0 0 0 0 B 0 71 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 mpi_comm_mct
S 2067 14 5 0 0 0 1 2063 13343 0 400000 A 0 0 0 0 B 0 71 0 0 0 0 0 480 3 0 0 0 0 0 0 0 0 0 0 0 0 71 0 624 0 0 0 0 cpl_attr_init cpl_attr_init 
F 2067 3 2064 2065 2066
A 166 2 0 0 0 6 1009 0 0 0 166 0 0 0 0 0 0 0 0 0 0 0
A 168 2 0 0 0 6 1010 0 0 0 168 0 0 0 0 0 0 0 0 0 0 0
A 339 2 0 0 0 6 1404 0 0 0 339 0 0 0 0 0 0 0 0 0 0 0
A 343 2 0 0 0 6 1405 0 0 0 343 0 0 0 0 0 0 0 0 0 0 0
A 345 2 0 0 0 6 1406 0 0 0 345 0 0 0 0 0 0 0 0 0 0 0
A 361 2 0 0 0 337 1413 0 0 0 361 0 0 0 0 0 0 0 0 0 0 0
A 362 2 0 0 0 337 1414 0 0 0 362 0 0 0 0 0 0 0 0 0 0 0
A 363 2 0 0 0 337 1415 0 0 0 363 0 0 0 0 0 0 0 0 0 0 0
A 364 2 0 0 0 337 1416 0 0 0 364 0 0 0 0 0 0 0 0 0 0 0
A 365 2 0 0 0 337 1417 0 0 0 365 0 0 0 0 0 0 0 0 0 0 0
A 366 2 0 0 0 337 1418 0 0 0 366 0 0 0 0 0 0 0 0 0 0 0
A 381 2 0 0 0 7 1409 0 0 0 381 0 0 0 0 0 0 0 0 0 0 0
A 382 2 0 0 0 7 1410 0 0 0 382 0 0 0 0 0 0 0 0 0 0 0
A 418 1 0 9 0 329 1448 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 426 1 0 9 0 339 1451 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 726 2 0 0 0 6 2030 0 0 0 726 0 0 0 0 0 0 0 0 0 0 0
Z
J 170 1 1
V 418 329 7 0
R 0 332 0 0
A 0 6 0 0 1 3 1
A 0 6 0 0 1 339 1
A 0 6 0 0 1 166 1
A 0 6 0 0 1 168 1
A 0 6 0 0 1 343 1
A 0 6 0 0 1 345 0
J 174 1 1
V 426 339 7 0
R 0 342 0 0
A 0 337 0 0 1 361 1
A 0 337 0 0 1 362 1
A 0 337 0 0 1 363 1
A 0 337 0 0 1 364 1
A 0 337 0 0 1 365 1
A 0 337 0 0 1 366 0
Z
