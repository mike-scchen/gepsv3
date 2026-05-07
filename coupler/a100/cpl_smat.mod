V34 :0x34 cpl_smat
12 cpl_smat.F90 S624 0
10/27/2025  15:42:25
use cpl_rank public 0 direct
use m_zeit private
use m_list private
use m_sparsematrixplus private
use m_sparsematrix private
enduse
D 341 26 1415 1840 1414 7
D 488 23 6 1 10 578 0 0 0 0 0
 10 577 11 10 577 578
D 491 23 6 1 11 578 0 0 0 0 0
 0 578 11 11 578 578
D 496 20 168
D 498 23 496 1 10 578 0 0 0 0 0
 10 577 11 10 577 578
D 501 23 496 1 11 578 0 0 0 0 0
 0 578 11 11 578 578
D 886 20 168
D 894 26 2383 11264 2382 7
S 624 24 0 0 0 10 1 0 5012 10005 8000 A 0 0 0 0 B 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 cpl_smat
S 626 23 0 0 0 10 1414 624 5036 4 0 A 0 0 0 0 B 400000 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 sparsematrix
S 627 19 0 0 0 10 1 624 5049 4 0 A 0 0 0 0 B 400000 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 262 1 0 0 0 0 0 624 0 0 0 0 sparsematrix_init
O 627 1 1530
S 629 19 0 0 0 10 1 624 5072 4 0 A 0 0 0 0 B 400000 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 260 1 0 0 0 0 0 624 0 0 0 0 sparsematrix_importgrowind
O 629 1 1609
S 631 19 0 0 0 10 1 624 5122 4 0 A 0 0 0 0 B 400000 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 258 1 0 0 0 0 0 624 0 0 0 0 sparsematrix_importgcolind
O 631 1 1615
S 633 19 0 0 0 10 1 624 5175 4 0 A 0 0 0 0 B 400000 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 256 2 0 0 0 0 0 624 0 0 0 0 sparsematrix_importmatelts
O 633 2 1639 1633
S 636 23 0 0 0 10 2382 624 5242 4 0 A 0 0 0 0 B 400000 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 sparsematrixplus
S 637 19 0 0 0 10 1 624 5259 4 0 A 0 0 0 0 B 400000 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 363 2 0 0 0 0 0 624 0 0 0 0 sparsematrixplus_init
O 637 2 2427 2417
S 638 19 0 0 0 10 1 624 5281 4 0 A 0 0 0 0 B 400000 9 0 0 0 0 0 0 0 0 0 0 0 0 0 0 360 1 0 0 0 0 0 624 0 0 0 0 sparsematrixplus_clean
O 638 1 2434
S 640 23 0 0 0 10 2402 624 5310 4 0 A 0 0 0 0 B 400000 10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 xonly
R 844 26 24 m_list =
S 1012 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 1013 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 1408 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
R 1414 25 3 m_sparsematrix sparsematrix
R 1415 5 4 m_sparsematrix nrows sparsematrix
R 1416 5 5 m_sparsematrix ncols sparsematrix
R 1417 5 6 m_sparsematrix data sparsematrix
R 1418 5 7 m_sparsematrix vecinit sparsematrix
R 1420 5 9 m_sparsematrix row_s sparsematrix
R 1421 5 10 m_sparsematrix row_s$sd sparsematrix
R 1422 5 11 m_sparsematrix row_s$p sparsematrix
R 1423 5 12 m_sparsematrix row_s$o sparsematrix
R 1425 5 14 m_sparsematrix row_e sparsematrix
R 1427 5 16 m_sparsematrix row_e$sd sparsematrix
R 1428 5 17 m_sparsematrix row_e$p sparsematrix
R 1429 5 18 m_sparsematrix row_e$o sparsematrix
R 1433 5 22 m_sparsematrix tcol sparsematrix
R 1434 5 23 m_sparsematrix tcol$sd sparsematrix
R 1435 5 24 m_sparsematrix tcol$p sparsematrix
R 1436 5 25 m_sparsematrix tcol$o sparsematrix
R 1440 5 29 m_sparsematrix twgt sparsematrix
R 1441 5 30 m_sparsematrix twgt$sd sparsematrix
R 1442 5 31 m_sparsematrix twgt$p sparsematrix
R 1443 5 32 m_sparsematrix twgt$o sparsematrix
R 1445 5 34 m_sparsematrix row_max sparsematrix
R 1446 5 35 m_sparsematrix row_min sparsematrix
R 1447 5 36 m_sparsematrix tbl_end sparsematrix
R 1530 14 119 m_sparsematrix init_
R 1609 14 198 m_sparsematrix importglobalrowindices_
R 1615 14 204 m_sparsematrix importglobalcolumnindices_
R 1633 14 222 m_sparsematrix importmatrixelementssp_
R 1639 14 228 m_sparsematrix importmatrixelementsdp_
S 1740 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 16 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 1741 3 0 0 0 6 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 1744 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 1745 3 0 0 0 7 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 7
S 1747 3 0 0 0 886 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10132 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 8 5b 4d 57 54 49 4d 45 5d
S 1748 3 0 0 0 886 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10141 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 8 5b 58 57 54 49 4d 45 5d
S 1749 3 0 0 0 886 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10150 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 8 5b 50 55 54 49 4d 45 5d
S 1750 3 0 0 0 886 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10159 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 8 5b 50 53 54 49 4d 45 5d
S 1751 3 0 0 0 886 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10168 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 8 5b 43 55 54 49 4d 45 5d
S 1752 3 0 0 0 886 1 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 10177 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 20 8 5b 43 53 54 49 4d 45 5d
R 1782 7 21 m_zeit masks$ac
R 1785 7 24 m_zeit header$ac
R 2382 25 4 m_sparsematrixplus sparsematrixplus
R 2383 5 5 m_sparsematrixplus strategy sparsematrixplus
R 2384 5 6 m_sparsematrixplus xprimelength sparsematrixplus
R 2385 5 7 m_sparsematrixplus xtoxprime sparsematrixplus
R 2386 5 8 m_sparsematrixplus yprimelength sparsematrixplus
R 2387 5 9 m_sparsematrixplus yprimetoy sparsematrixplus
R 2388 5 10 m_sparsematrixplus matrix sparsematrixplus
R 2389 5 11 m_sparsematrixplus tag sparsematrixplus
R 2402 16 24 m_sparsematrixplus xonly
R 2417 14 39 m_sparsematrixplus initfromroot_
R 2427 14 49 m_sparsematrixplus initdistributed_
R 2434 14 56 m_sparsematrixplus clean_
S 2457 6 4 0 0 341 2458 624 13365 4 8 A 0 0 0 0 B 0 14 0 0 0 0 0 0 0 0 0 0 2470 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfs2rsm_smat
S 2458 6 4 0 0 341 2459 624 13378 4 8 A 0 0 0 0 B 0 15 0 0 0 1840 0 0 0 0 0 0 2470 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfs2gocn_smat
S 2459 6 4 0 0 341 2460 624 13392 4 8 A 0 0 0 0 B 0 16 0 0 0 3680 0 0 0 0 0 0 2470 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocn2gfs_smat
S 2460 6 4 0 0 341 2461 624 13406 4 8 A 0 0 0 0 B 0 17 0 0 0 5520 0 0 0 0 0 0 2470 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rsm2rocn_smat
S 2461 6 4 0 0 341 2462 624 13420 4 8 A 0 0 0 0 B 0 18 0 0 0 7360 0 0 0 0 0 0 2470 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rocn2rsm_smat
S 2462 6 4 0 0 341 2463 624 13434 4 8 A 0 0 0 0 B 0 19 0 0 0 9200 0 0 0 0 0 0 2470 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocn2rocn_smat
S 2463 6 4 0 0 894 2464 624 13449 4 8 A 0 0 0 0 B 0 21 0 0 0 11040 0 0 0 0 0 0 2470 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfs2rsm_smatp
S 2464 6 4 0 0 894 2465 624 13463 4 8 A 0 0 0 0 B 0 22 0 0 0 22304 0 0 0 0 0 0 2470 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfs2gocn_smatp
S 2465 6 4 0 0 894 2466 624 13478 4 8 A 0 0 0 0 B 0 23 0 0 0 33568 0 0 0 0 0 0 2470 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocn2gfs_smatp
S 2466 6 4 0 0 894 2467 624 13493 4 8 A 0 0 0 0 B 0 24 0 0 0 44832 0 0 0 0 0 0 2470 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rsm2rocn_smatp
S 2467 6 4 0 0 894 2468 624 13508 4 8 A 0 0 0 0 B 0 25 0 0 0 56096 0 0 0 0 0 0 2470 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rocn2rsm_smatp
S 2468 6 4 0 0 894 1 624 13523 4 8 A 0 0 0 0 B 0 26 0 0 0 67360 0 0 0 0 0 0 2470 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocn2rocn_smatp
S 2469 19 0 0 0 6 1 624 5067 4 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 368 3 0 0 0 0 0 624 0 0 0 0 init
O 2469 3 2427 2417 1530
S 2470 11 0 0 0 10 2456 624 13539 40800000 805000 A 0 0 0 0 B 0 28 0 0 0 78624 0 0 2457 2468 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _cpl_smat$0
S 2471 23 5 0 0 0 2478 624 13551 0 0 A 0 0 0 0 B 0 29 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 cpl_smat_init
S 2472 1 3 1 0 6 1 2471 13565 4 3000 A 0 0 0 0 B 0 29 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 myrank
S 2473 1 3 1 0 6 1 2471 13572 4 3000 A 0 0 0 0 B 0 29 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rank_root
S 2474 1 3 1 0 6 1 2471 13582 4 3000 A 0 0 0 0 B 0 29 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 id_src
S 2475 1 3 1 0 6 1 2471 13589 4 3000 A 0 0 0 0 B 0 29 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 id_dst
S 2476 1 3 1 0 6 1 2471 13596 4 3000 A 0 0 0 0 B 0 29 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 mpi_comm_mct
S 2477 1 3 1 0 30 1 2471 13609 4 43000 A 0 0 0 0 B 0 29 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 scrip_file
S 2478 14 5 0 0 0 1 2471 13551 0 400000 A 0 0 0 0 B 0 29 0 0 0 0 0 676 6 0 0 0 0 0 0 0 0 0 0 0 0 29 0 624 0 0 0 0 cpl_smat_init cpl_smat_init 
F 2478 6 2472 2473 2474 2475 2476 2477
S 2479 23 5 0 0 0 2483 624 13620 0 0 A 0 0 0 0 B 0 89 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 read_scrip_file
S 2480 1 3 1 0 30 1 2479 13609 4 43000 A 0 0 0 0 B 0 89 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 scrip_file
S 2481 1 3 2 0 341 1 2479 13636 4 3000 A 0 0 0 0 B 0 89 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 remap_smat
S 2482 1 3 1 0 6 1 2479 12902 4 3000 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 tag
S 2483 14 5 0 0 0 1 2479 13620 0 400000 A 0 0 0 0 B 0 89 0 0 0 0 0 683 3 0 0 0 0 0 0 0 0 0 0 0 0 89 0 624 0 0 0 0 read_scrip_file read_scrip_file 
F 2483 3 2480 2481 2482
A 166 2 0 0 0 6 1012 0 0 0 166 0 0 0 0 0 0 0 0 0 0 0
A 168 2 0 0 0 6 1013 0 0 0 168 0 0 0 0 0 0 0 0 0 0 0
A 343 2 0 0 0 6 1408 0 0 0 343 0 0 0 0 0 0 0 0 0 0 0
A 539 2 0 0 0 6 1740 0 0 0 539 0 0 0 0 0 0 0 0 0 0 0
A 541 2 0 0 0 6 1741 0 0 0 541 0 0 0 0 0 0 0 0 0 0 0
A 557 2 0 0 0 496 1747 0 0 0 557 0 0 0 0 0 0 0 0 0 0 0
A 558 2 0 0 0 496 1748 0 0 0 558 0 0 0 0 0 0 0 0 0 0 0
A 559 2 0 0 0 496 1749 0 0 0 559 0 0 0 0 0 0 0 0 0 0 0
A 560 2 0 0 0 496 1750 0 0 0 560 0 0 0 0 0 0 0 0 0 0 0
A 561 2 0 0 0 496 1751 0 0 0 561 0 0 0 0 0 0 0 0 0 0 0
A 562 2 0 0 0 496 1752 0 0 0 562 0 0 0 0 0 0 0 0 0 0 0
A 577 2 0 0 0 7 1744 0 0 0 577 0 0 0 0 0 0 0 0 0 0 0
A 578 2 0 0 0 7 1745 0 0 0 578 0 0 0 0 0 0 0 0 0 0 0
A 614 1 0 9 0 488 1782 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 622 1 0 9 0 498 1785 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
J 170 1 1
V 614 488 7 0
R 0 491 0 0
A 0 6 0 0 1 3 1
A 0 6 0 0 1 343 1
A 0 6 0 0 1 166 1
A 0 6 0 0 1 168 1
A 0 6 0 0 1 539 1
A 0 6 0 0 1 541 0
J 174 1 1
V 622 498 7 0
R 0 501 0 0
A 0 496 0 0 1 557 1
A 0 496 0 0 1 558 1
A 0 496 0 0 1 559 1
A 0 496 0 0 1 560 1
A 0 496 0 0 1 561 1
A 0 496 0 0 1 562 0
Z
