# timcom_glb

## Quick Steup

This guide walks you through building the GEPS library and running tests on both ARM and X86 platforms.

### 1. Build the GEPS Library

```sh
./build.sh <platform> # platform = fx1000 | a100
```

### 2. Run OMIP

OMIP validates ocean model accuracy and performance.

```sh
./build.sh <platform> OMIP # platform = fx1000 | a100
```

After building, navigate to the jobs directory and submit the job:

```sh
cd jobs
./build_job.sh <platform>  # platform = fx1000 | a100
cd 20010100
pjsub submit_job.sh
```

### 3. Run unit tests on X86 (HGX-A100)

Unit tests ensure that the results between the CPU (X86) and GPU (A100) are consistent.

```sh
./build.sh a100 OMIP
cd jobs
./build_job.sh a100-qa
cd 20010100
pjsub submit_job.sh
```
