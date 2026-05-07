#/bin/bash

FCT_MODEL=${kshpath}/../build_a100/src/timcom_omip

cat > ${MDIR}/submit_job.sh << EOF
#!/bin/bash
#PJM -L vnode=1
#PJM -L vnode-core=32
#PJM -L gpu=8
#PJM -L rscunit=rscunit_pg01
#PJM -L rscgrp=gpu-rd-large
#PJM -x PJM_CACHE_MODE=4
#PJM -L elapse=03:00:00
#PJM -L node-mem=unlimited
#PJM --no-stging
#PJM --mpi "proc=${MPI}"
#PJM -j
#PJM -g sum
#PJM -N op_${dtg}
#PJM -o %j.log
#PJM -e %j.err

. /usr/share/Modules/init/bash
module use /package/x86_64/nvidia/hpc_sdk/modulefiles
module load nvhpc-hpcx-cuda12/24.11
module unuse /package/x86_64/nvidia/hpc_sdk/modulefiles

export MPI=${MPI}
export NVCOMPILER_ACC_CUDA_NOCOPY=1
export NVCOMPILER_ACC_CUDA_MEMALLOCASYNC="1"
export NVCOMPILER_ACC_CUDA_MEMALLOCASYNC_POOLSIZE="64G"
export NVCOMPILER_ACC_USE_GRAPH="1"

cd ${MDIR}

FCT_MODEL=${FCT_MODEL}
/usr/bin/time -p mpiexec -n \${MPI} \${FCT_MODEL} 000

if [ \$? != 0 ] ; then
  echo "error occured: fct model fail !!"
fi

echo "Ending at: " \$(date)
EOF
