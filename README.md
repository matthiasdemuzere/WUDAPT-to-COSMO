# WUDAPT-to-COSMO

Set of tools to use Local Climate Zone (LCZ)-based urban canopy parameters in DWD's COSMO-CLM NWP and regional climate model.

## Citaton
Varentsov, M., Samsonov, T., Demuzere, M., (2020). Impact of urban canopy parameters on a megacity’s modelled thermal environment. Atmosphere 11(12), 1349; [https://doi.org/10.3390/atmos11121349](https://www.mdpi.com/2073-4433/11/12/1349).

## Context
TERRA_URB is the urban canopy parameterization embedded in TERRA-ML, the land surface model in COSMO-CLM. By default it uses impervious surface area information from the [Copernicus Land Monitoring Service](https://land.copernicus.eu/pan-european/high-resolution-layers/imperviousness) (for Europe) / [National Geophysical Data Center](https://databasin.org/datasets/016d2235a5ed43ad83ceeed6c408d149) (global) and anthropogenic heat flux information from [Flanner et al. (2010)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008gl036465). All other geometrical, thermal and radiative urban canopy parameters are spatially invariant, and set to the bulk values provided in Table 1 of [Wouters et al., 2016](https://gmd.copernicus.org/articles/9/3027/2016/)).

The set of tools provided in this repo allow one to introduce LCZ-based urban canopy parameters, compiled from [Stewart and Oke (2012)](http://10.1175/BAMS-D-11-00019.1) and [Stewart et al. (2014)](http://10.1002/joc.3746).

This work is an outcome of AEVUS I and II, the COSMO Priority Tasks on "Analysis and evaluation of the TERRA_URB scheme". More info [here](http://www.cosmo-model.org/content/tasks/priorityTasks/default.htm) (project pages only accessible to COSMO members). Preliminary test results of LCZ parameters in COSMO-CLM are also described in Brousse et al. ([2019](https://doi.org/10.1016/j.uclim.2018.12.004), [2020](https://onlinelibrary.wiley.com/doi/abs/10.1002/joc.6477)) and [Van de Walle et al. (2021)](http://doi.org/10.1007/s00704-021-03733-7).  



## Requirements
* Be a member of the [COSMO-CLM community](https://wiki.coast.hzg.de/clmcom/), in order to be able to access [EXTPAR](https://wiki.coast.hzg.de/clmcom/external-data-98599196.html).
* Have your domain file available from EXTPAR (netcdf file)
* an LCZ map covering the same region of interest. Sources for existing LCZ maps are listed [here](https://www.wudapt.org/lcz-maps/).


## Virtual environment

It is advised to use a python virtual environment, eg.:

```
> python3.9 -m venv ~/w2c_venv
> . ~/w2c_venv/bin/activate
> pip install -r requirements.txt
```

The `requirements.txt` can be generated using `pipreqs`: 
```
pipreqs --ignore=terra/ .
```


### Execute

The run code is currently configured for the Moscow case, as developed in Varentsov et al.
CLM and LCZ input data used in this study is provided under `data/`.

```
python run.py
```

