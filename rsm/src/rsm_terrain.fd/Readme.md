# RSM Variant of the `mtnlm7_slm30g` Program

## Changelog

* Remove Gaussian grid generation as it is inappropriate for regional domain.
* Fix computation error due to the inappropriate assumption of periodic boundary for regional domain. For example, orographic asymmetry (OA) is affected.
* Implement new map projection functionality to generate grid coordinates for regional domain. Mercator and polar stereographic are supported as in RSM.
* Implement new data selection functionality to perform computations only within regional domain.
* Implement new I/O functionality to read map projection parameters directly from namelist.
* Cleanup unused code.
* Drop dependency on NCEPLIBS.
* Simplify Makefile.
* Simplify run scripts.

## Compilation

```bash
make -f Makefile_slm30g_rsm clean
make -f Makefile_slm30g_rsm
```

## Usage

1. Modify the map projection parameters in the `rsm-domain.nml` namelist as needed.
2. Execute directly in shell:

```bash
./mtnlm7_slm30g.exe <<< "IGRD JGRD"
```

where `IGRD`: Domain zonal dimension; `JGRD`: Domain meridional dimension.

3. Alternatively, modify the `run_mtnlm7_slm30g.sh` shell script and submit it to the job scheduler for execution.
4. The following output files should be generated:
	* `rsm-lat.bin`: Latitude coordinates for each grid point.
	* `rsm-lon.bin`: Longitude coordinates for each grid point.
	* `rsm-oro.bin`: Topography.
	* `rsm-oro-smoothed.bin`: Smoothed topography.
	* `rsm-lsm.bin`: Land sea mask.
	* `rsm-hprime.bin`: 14 orographic variables, referred to as <i>h<sup>′</sup></i>.
