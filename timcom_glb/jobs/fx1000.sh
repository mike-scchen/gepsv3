#/bin/bash

FCT_MODEL=${kshpath}/../build_fx1000/src/timcom_omip

cat > ${MDIR}/submit_job.sh << EOF
#!/bin/bash
#PJM -L "node=12:noncont"
#PJM -L rscgrp=small
#PJM -x PJM_CACHE_MODE=4
#PJM -L elapse=3:00:00
#PJM -L node-mem=unlimited
#PJM --no-stging
#PJM --mpi "proc=${MPI}"
#PJM -j
#PJM -g sum
#PJM -N op_${dtg}
#PJM -o %j.log
#PJM -e %j.err

OMP=1
source /users/xa09/sample/setup_mpi+omp.fx1000 \$OMP

export MPI=${MPI}

cd ${MDIR}

FCT_MODEL=${FCT_MODEL}
/usr/bin/time -p mpiexec -n \${MPI} \${FCT_MODEL} 000

if [ \$? != 0 ] ; then
  echo "error occured: fct model fail !!"
fi

echo "Ending at: " \$(date)
EOF
