V34 :0x34 cpl_gsmap
13 cpl_gsmap.F90 S624 0
10/27/2025  15:42:24
use cpl_rank public 0 direct
use m_globalsegmap private
enduse
D 58 26 646 472 643 7
D 185 23 6 1 11 172 0 0 1 0 0
 0 171 11 11 172 172
D 188 23 6 1 11 172 0 0 1 0 0
 0 171 11 11 172 172
S 624 24 0 0 0 10 1 0 5012 10005 8000 A 0 0 0 0 B 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 cpl_gsmap
S 626 23 0 0 0 10 643 624 5037 4 0 A 0 0 0 0 B 400000 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 globalsegmap
S 627 19 0 0 0 10 1 624 5050 4 0 A 0 0 0 0 B 400000 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 66 6 0 0 0 0 0 624 0 0 0 0 globalsegmap_init
O 627 6 805 797 787 768 747 727
S 629 19 0 0 0 10 1 624 5073 4 0 A 0 0 0 0 B 400000 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 59 1 0 0 0 0 0 624 0 0 0 0 globalsegmap_lsize
O 629 1 853
S 631 19 0 0 0 10 1 624 5098 4 0 A 0 0 0 0 B 400000 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 57 1 0 0 0 0 0 624 0 0 0 0 globalsegmap_clean
O 631 1 813
R 643 25 1 m_globalsegmap globalsegmap
R 646 5 4 m_globalsegmap comp_id globalsegmap
R 647 5 5 m_globalsegmap gsize globalsegmap
R 652 5 10 m_globalsegmap ngseg globalsegmap
R 666 5 24 m_globalsegmap start globalsegmap
R 667 5 25 m_globalsegmap start$sd globalsegmap
R 668 5 26 m_globalsegmap start$p globalsegmap
R 669 5 27 m_globalsegmap start$o globalsegmap
R 672 5 30 m_globalsegmap length globalsegmap
R 673 5 31 m_globalsegmap length$sd globalsegmap
R 674 5 32 m_globalsegmap length$p globalsegmap
R 675 5 33 m_globalsegmap length$o globalsegmap
R 678 5 36 m_globalsegmap pe_loc globalsegmap
R 679 5 37 m_globalsegmap pe_loc$sd globalsegmap
R 680 5 38 m_globalsegmap pe_loc$p globalsegmap
R 681 5 39 m_globalsegmap pe_loc$o globalsegmap
R 727 14 85 m_globalsegmap initd_
R 747 14 105 m_globalsegmap initr_
R 768 14 126 m_globalsegmap initp_
R 787 14 145 m_globalsegmap initp1_
R 797 14 155 m_globalsegmap initp0_
R 805 14 163 m_globalsegmap init_index_
R 813 14 171 m_globalsegmap clean_
R 853 14 211 m_globalsegmap lsize_
S 962 6 4 0 0 58 963 624 6318 4 8 A 0 0 0 0 B 0 9 0 0 0 0 0 0 0 0 0 0 974 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfsgsmap
S 963 6 4 0 0 58 964 624 6327 4 8 A 0 0 0 0 B 0 10 0 0 0 472 0 0 0 0 0 0 974 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocn_in_gfsgsmap
S 964 6 4 0 0 58 965 624 6344 4 8 A 0 0 0 0 B 0 11 0 0 0 944 0 0 0 0 0 0 974 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rsm_in_gfsgsmap
S 965 6 4 0 0 58 966 624 6360 4 8 A 0 0 0 0 B 0 13 0 0 0 1416 0 0 0 0 0 0 974 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rsmgsmap
S 966 6 4 0 0 58 967 624 6369 4 8 A 0 0 0 0 B 0 14 0 0 0 1888 0 0 0 0 0 0 974 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rocn_in_rsmgsmap
S 967 6 4 0 0 58 968 624 6386 4 8 A 0 0 0 0 B 0 15 0 0 0 2360 0 0 0 0 0 0 974 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfs_in_rsmgsmap
S 968 6 4 0 0 58 969 624 6402 4 8 A 0 0 0 0 B 0 17 0 0 0 2832 0 0 0 0 0 0 974 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocngsmap
S 969 6 4 0 0 58 970 624 6412 4 8 A 0 0 0 0 B 0 18 0 0 0 3304 0 0 0 0 0 0 974 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gfs_in_gocngsmap
S 970 6 4 0 0 58 971 624 6429 4 8 A 0 0 0 0 B 0 19 0 0 0 3776 0 0 0 0 0 0 974 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rocn_in_gocngsmap
S 971 6 4 0 0 58 972 624 6447 4 8 A 0 0 0 0 B 0 21 0 0 0 4248 0 0 0 0 0 0 974 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rocngsmap
S 972 6 4 0 0 58 973 624 6457 4 8 A 0 0 0 0 B 0 22 0 0 0 4720 0 0 0 0 0 0 974 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 rsm_in_rocngsmap
S 973 6 4 0 0 58 1 624 6474 4 8 A 0 0 0 0 B 0 23 0 0 0 5192 0 0 0 0 0 0 974 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 gocn_in_rocngsmap
S 974 11 0 0 0 10 961 624 6492 40800000 805000 A 0 0 0 0 B 0 25 0 0 0 5664 0 0 962 973 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _cpl_gsmap$0
S 975 23 5 0 0 0 984 624 6505 0 0 A 0 0 0 0 B 0 26 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 cpl_gsmap_init
S 976 1 3 1 0 6 1 975 6520 4 3000 A 0 0 0 0 B 0 26 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 myrank
S 977 6 3 1 0 6 1 975 6527 800004 3000 A 0 0 0 0 B 0 26 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 seg_size
S 978 7 3 1 0 185 1 975 6536 800204 3000 A 0 0 0 0 B 0 26 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 seg_strt
S 979 7 3 1 0 188 1 975 6545 800204 3000 A 0 0 0 0 B 0 26 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 seg_leng
S 980 1 3 1 0 6 1 975 6554 4 3000 A 0 0 0 0 B 0 26 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 rank_root
S 981 1 3 1 0 6 1 975 6564 4 3000 A 0 0 0 0 B 0 26 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 mct_comm_world
S 982 1 3 1 0 6 1 975 6579 4 3000 A 0 0 0 0 B 0 26 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 id_comp
S 983 1 3 0 0 30 1 975 6587 4 43000 A 0 0 0 0 B 0 26 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 msg
S 984 14 5 0 0 0 1 975 6505 200 400000 A 0 0 0 0 B 0 26 0 0 0 0 0 118 8 0 0 0 0 0 0 0 0 0 0 0 0 26 0 624 0 0 0 0 cpl_gsmap_init cpl_gsmap_init 
F 984 8 976 977 978 979 980 981 982 983
S 985 6 1 0 0 7 1 975 6591 40800006 3000 A 0 0 0 0 B 0 30 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 z_e_171
A 170 1 0 0 0 6 977 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 171 7 0 0 0 7 170 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
A 172 1 0 0 0 7 985 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
Z
Z
