kshpath=`pwd `
. ${kshpath}/00_setDate.sh

 
 Caldtg=/data/common/gfs/scripts/Caldtg.ksh
 cd ${kshpath}
#######################################################################
# 3. Set date and update hours
#######################################################################

  date_ini=${pred_str}
  updatehr=${gap}
  date_end=${pred_end}
  echo date_end=${date_end}
  declare -a Case_a=( "000" "001" "002" "003" "004" "005" "006" "007" "008"  \
	  "009" "010" "011" "012" "013" "014" "015" "016" "017" "018" "019" \
	  "020" "021" "022" "023" "024" "025" "026" "027" "028" "029" "030" \
	  "031" "032")


tmpdtg=${date_ini}

  while [ ${tmpdtg} -le ${date_end} ]
  do
    dtg_d2=`echo ${tmpdtg} | cut -c3-10`
    for mem in $(seq 0 1 $member) ;
    do
      if [ -f ${kshpath}/${dtg_d2}_${Case_a[$mem]}.ufs/rsm/exp_run.sh ]; then
         echo " ${dtg_d2} finished "
      else
         echo "build_job.sh ${dtg_d2}"
         echo ${tmpdtg} > ${kshpath}/crdate
         if [ ${em_det} -eq 0 ]; then
#---------------build det---------------------------------------------------------
         export suffix="${Case_a[$em_det]}"
/usr/bin/newgrp sum << eof
           bash build_job.sh
eof
 	 else
#---------------build em---------------------------------------------------------
              export suffix="${Case_a[$mem]}"
/usr/bin/newgrp sum << eof
              cp ${kshpath}/config_gfs_fx1000_em.sh_ ${kshpath}/config_gfs_fx1000_em.sh
              sed -i 's/MMM/${Case_a[$mem]}/g' ${kshpath}/config_gfs_fx1000_em.sh
	      bash build_job.sh
eof
	 fi

#---------------Check Ocean IC---------------------------------------------------------
	 if [ -f /data/common/gfs/GEPSv3_lib/data/HYCOM_ic/TIMCOM_g1536_${tmpdtg}.nc ]; then
              echo "check ic ok :"${tmpdtg} >> ${kshpath}/check_ic.txt
           else
              chk=`${Caldtg} ${tmpdtg} -24`
              chk2=${tmpdtg}
              export chkic=`echo ${chk} |cut -c3-10`
              export chkic2=`echo ${chk2} |cut -c3-10`
#/usr/bin/newgrp sum << eof
#              /package/x86_64/ncl-6.6.2/bin/ncl check_ic.ncl
#eof
	      echo "check ic false :"${tmpdtg} >> ${kshpath}/check_ic.txt
         fi

	 if [ ${em_det} -eq 0 ]; then
#---------------submit det---------------------------------------------------------
           echo "Auto-submission disabled."
	   echo "submit the job by: pjsub (date)_000.ufs/submit_job.sh"
           #/usr/bin/pjsub ${dtg_d2}_000.ufs/submit_job.sh
         else
#---------------submit em---------------------------------------------------------
           echo "Auto-submission disabled."
	   echo "submit the job by: pjsub (date)_(member).ufs/submit_job.sh"
           #/usr/bin/pjsub ${dtg_d2}_${Case_a[$mem]}.ufs/submit_job.sh
         fi
      fi
    done
    tmpdtg=`${Caldtg} ${tmpdtg} ${updatehr}`
  done


