# LCZ-to-CCLM

Set of tools to use Local Climate Zone (LCZ)-based urban canopy parameters in DWD's COSMO-CLM NWP and regional climate model.

## Context
TERRA_URB, the urban canopy parameterization in COSMO-CLM ([Wouters et al., 2016](https://gmd.copernicus.org/articles/9/3027/2016/)), by default uses impervious surface area information from the [Copernicus Land Monitoring Service](https://land.copernicus.eu/pan-european/high-resolution-layers/imperviousness) (for Europe) / [National Geophysical Data Center](https://databasin.org/datasets/016d2235a5ed43ad83ceeed6c408d149) (global) and anthropogenic heat flux information from [Flanner et al. (2010)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008gl036465). All other geometrical, thermal and radiative urban canopy parameters are spatially invariant, and set to the values provided in Table 1 of [Wouters et al., 2016](https://gmd.copernicus.org/articles/9/3027/2016/)).

The set of tools provided in this repo allow one to introduce LCZ-based urban canopy parameters, compiled from [Stewart and Oke (2012)](http://10.1175/BAMS-D-11-00019.1) and [Stewart et al. (2014)](http://10.1002/joc.3746).



## Requirements
* Be a member of the [COSMO-CLM community](https://wiki.coast.hzg.de/clmcom/), in order to be able to access [EXTPAR](https://wiki.coast.hzg.de/clmcom/external-data-98599196.html).
* have your domain file available from EXTPAR (netcdf file)
* an LCZ map covering the same region of interest. If not available, please contact me.


## Procedure (more info to come)

`from utils import *`

### Step 1: Assign UCP values to LCZ map and convert to COSMO Grid
`lcz_to_cosmo(ucpFile, clmFile, lczFile, bandNr, ucpVersion, nrLcz=17,
              interpMethod='linear', aggregation=True, aggregationScale=2,
              isaWeight=True, saiWeight=False, fileNameExt='')`

### Step 2: Address the double counting issue.
`removeDoubleCounting(clmFileNew,gcFile,removeUrban=False)`



## Python version, required packages
This code was developed and tested using python 3.6.10.

* pandas
* numpy
* xarray
* scipy


