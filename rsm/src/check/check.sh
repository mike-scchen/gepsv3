cd /nwpr/gfs/xb157/CWBSUMc2/4cpl/op_work_4cpl_typ/rsm/wrk/dir_cmp/cmp_rsm_4cpl/
ff=/nwpr/gfs/xb157/CWBSUMc2/4cpl/op_work_4cpl_typ/rsm/
for fl in `ls *.F`
do
	echo "${fl}"
	diff ${fl} ~xb119/tcogit/cwbsum/rsm/wrk/dir_cmp/cmp_rsm/${fl} > ${ff}src/check/${fl}_diff
done

for fl in `ls *.f90`
do
        echo "${fl}"
        diff ${fl} ~xb119/tcogit/cwbsum/rsm/wrk/dir_cmp/cmp_rsm/${fl} > ${ff}src/check/${fl}_diff
done

