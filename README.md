# GEPSv3op3.3 (plz git clone at dm16)

  ## 1.get
  ```sh
  git config --global http.sslVerify false
  git config --global http."https://git.rdc.cwb/".sslVerify false
  git clone --recursive git@git.rdc.cwb:tco/GEPSv3.git
  ```
  TCO:4648ba08c82267c702b2fba056773fac564977c4
  TIMCOM:bec340aaab12008a8ceaea46762c5c0a45d93407
  RSM:a84a5f6e4c75c1fdb15e9fd1a505d6706548e49f

  ## 2.compile
  ### on a100 (GPU) - only (Extended weather)
  only glb 2cpl
  ```sh
    cd tco639l72/
    ./compile a100 to make build_a100/lib/libtcogfs.a file
      #CPL ( vi compile, TIMCOMCPL=ON ) !check
      #default is uncpl
      
    cd timcom_glb/
    ./build.sh a100 to make build_a100/src/libtimcom.a file

    cd driver
    ./build.sh a100 2cpl TCo383/TCo199 to make exec 
           "TCoTIMCOM_2cpl_a100_TCo383/TCo199"

    cd jobs
    vi 00_setDate.sh set
      export JCAP="383" #383 639
      export mach="a100" #a100 fx1000
      export struc="2cpl" #2cpl(a100) 4cpl(fx1000)

    vi config_gfs.sh
    ./99_main.sh (submit jobs)
  ```
  ### on fx1000 (CPU) 2cpl_CICE, 2cpl, 4cpl
  #### **cice 
  (Extended weather can skip this step)
  ```sh
      cd cice/
      ./run_cice_build.sh
  ```

  #### **TCo
  ```sh
      cd tco639l72/
  ```
  ##### use "build" - Suggested to use
  ```sh
      ./build.sh <machine> <option> ;
      ./build.sh fx1000      (only atm)
      ./build.sh fx1000 2cpl (glb cpl)
      ./build.sh fx1000 4cpl (all)
      ./build.sh fx1000 2cpl_CICE (glb cpl + cice)

  ```
  ##### use "cmake" - only (Extended weather)
  ```sh
      ./compile fx1000 
        #default is uncpl
        #2CPL (vi compile, -DTIMCOMCPL=ON)
        #4CPL (vi compile, -DUSE_RSM=ON -DCWBSUM=ON -DTIMCOMCPL=ON)
  ```
  #### **TIMCOM glb
  ```sh
      cd timcom_glb/
      vi build.sh -> USE_CLM=OFF (Extended weather), USE_CLM=ON (CLM+CICE)
      ./build.sh fx1000
  ```
  #### **Regional model
  ##### !! The following is only needed for 4cpl; skip for 2cpl(glb cpl) !!
  ```sh
      cd rsm/
      ./build.sh fx1000 5km

      cd rsm/src/timcom.fd
      make
  ```

  ##### check rsm/run/exe/librsm4cpl.a exist
  ## 3.get exec
  ```sh
  cd driver
  ./build.sh [machine] [frame] [ATM resolution 383/199 ]
  ./build.sh a100 2cpl TCo383  to make glb_cpl(GPU) exec
  ./build.sh fx1000 2cpl TCo383 to make glb_cpl(CPU) exec 
  ./build.sh fx1000 4cpl TCo383 to make GEPSv3 exec
  ./build.sh fx1000 2cpl_CICE TCo383
  ```
  ## 4.run
  ```sh
  cd jobs
  vi 00_setDate.sh set init time fcsthr
     export JCAP="383" #383 199
     export mach="fx1000" #a100 fx1000
     export struc="4cpl" #2cpl(a100) 2cpl_CICE、2cpl、4cpl(fx1000) 
  
  ./99_main.sh (submit jobs)
  ```
  ### run CLM (restart) 
  ```sh
  need rebuild timcom_glb
    cd timcom_glb/
    vi build.sh
    USE_CLM=ON
    continue 4.get exec and 5.run
  cd jobs 
  vi config_gfs*
     dorst=true
  ```
  #### do-restart
  ```sh
  check wrk/gfsctl not begin from 0 (from 0 = no-restart)
        maybe set from 24
  check wrk/timectl restart file output time > your restart time
        >24
  check wrk/namelist.run restart=true, set init_time(begin time) and duration
  ```




