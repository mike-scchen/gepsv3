kshpath=` pwd `

export check_pre=`date -d "-1 days" +%Y%m%d`00
export check=`echo ${check_pre} |cut -c3-10`
export op_date=`cat ${kshpath}/crdate |cut -c1-10`

#set run dtg time
#export pred_str=`date -d "-3 days" +%Y%m%d`00 #2023121000 #${hrname}
#export pred_end=`date -d "-3 days" +%Y%m%d`00 #2023121000
export pred_str=2025010100 #${hrname}
export pred_end=2025010100
#set interval
export gap=24           

export RESTRHR=0        #output initial time
export FCSTGAP=24        #output frequency (ocn)-> ../tco639l72/etc_4cpl/ocards (atm)
export FCSTHR=24       #forecast hour-> Must be consistent with ../tco639l72/etc_4cpl/gfsctl

export FRE_R=24        #restart output freq.

export ft_day=$(( ${FCSTHR}/${gap} ))
export ice_year=`cat ${kshpath}/crdate |cut -c1-4`
export ice_mon=`cat ${kshpath}/crdate |cut -c5-6`
export ice_day=`cat ${kshpath}/crdate |cut -c7-8`

#set det or em
export em_det=0  # 1=em, 0=det
if [ ${em_det} -eq 1 ]; then
  export member=2   #set how many member
else
  export member=0
fi

#set something
export JCAP="383" #199 383 639
export mach="fx1000" #a100 fx1000
export struc="2cpl_CICE" #2cpl(gpu), 2cpl / 4cpl / 2cpl_CICE(fx1000) 
export tai_v='_op'  #Regional Version
export FCT_MODEL=${kshpath}/../driver/TCoTIMCOM_${struc}_${mach}_TCo${JCAP}
