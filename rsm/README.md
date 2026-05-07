# This is the reame file for RSM/MSM model (git@git.rdc.cwb:yjchen/rsm.git).

  RSM is hydrostatic Regional Spectral Model [Juang and Kanamitsu (1994), Juang et al. (1997)] and MSM (Mesoscale Spectral Model; Juang (2000)) is the non-hydrostatic version.

  _written by Chen, Ying-Ju 2021.08.02_


## model directories and files

   README.md : This file
   build.sh  : link/build libraries, data, post-process exe
   fix/      : fix-data
   job/      : scripts for running model at HPC
   lib/      : libraries
   opt/      : default options
   run/      : model setting
   script/   : scripts for creating terrain, LSM data, etc.
   src/      : source code
   utl/      : utilities
   wrk/      : working directory (will be created after build.sh),
               including the compiled code, GFS base field, and model output


## How to build?

   There are 5 options:
   1. a default test case, type:
   `./build.sh [machine] test`

   2. To get default RSM 5km setting:
   `./build.sh [machine] 5km`

   3. To get default RSM 12km setting:
   `./build.sh [machine] 12km`

   4. To setup your own setting:
   `./build.sh [machine] user`

   5. build without setting any experiment:
   `./build.sh [machine] no`
    In this case, only the libraries, utilies and fix data are linked.

   The supported machines are fx10, fx100, pcc.


## How to run?

   If the model is built and ready to be executed,
   go to job/fujitsu (fx10/fx100) and type
   `pjsub run_rsm.sh`
   or
   go to job/pcc and type
   `qsub run_rsm.sh`


## Post-process

   To convert binary output to grib1 files after the 
   model finishes the integration,
   go to job/fujitsu (fx10/fx100) and type
   `pjsub run_post.sh`
   or
   go to job/pcc and type
   `qsub run_post.sh`

   The default output directory is at rsm/wrk/${SDATE}


## Note
   
   1. If you want to change the domain, clean the old executable files and terrain data first!
      ```
      cd run/
      ./clean
      ```

   2. Some changes in run/configure will required re-compilation. 
      To know what options need to re-compile, check run/compile to find the keyworkd, and then grep the key word in src/. 
      If there are macro-define (`#ifdef KEYWORD`) of the key words in src/, then whenever you change the related options in run/cofigure, you need to compile again to make the changes work!
      * For example:
        If you have done `./build.sh [at_MACHINE] [do_case]`, which means you have
        - [x] linked the basefield data from GFS
        - [x] linked the libraries
        - [x] built the post-process executable (utl/rpgbnawips.x)
        and you want to re-compile RSM/MSM,
        you can go to run/ `cd run` and type `./compile` to setup your new setting.

   3. To debug, it is easier to check the source code after compilation (renamed to \*.f90).
      For exmaple, if anything is wrong in the forecasting integration part, go check the source code in wrk/dir_cmp/cmp_rsm/\*.f90


  
