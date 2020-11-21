# WUDAPT-to-COSMO

Set of tools to use Local Climate Zone (LCZ)-based urban canopy parameters in DWD's COSMO-CLM NWP and regional climate model.

## Citaton
Varentsov, M., Samsonov, T., Demuzere, M., (under review). Impact of urban canopy parameters on a megacity’s modelled thermal environment. Atmosphere.

## Context
TERRA_URB is the urban canopy parameterization embedded in TERRA-ML, the land surface model in COSMO-CLM. By default it uses impervious surface area information from the [Copernicus Land Monitoring Service](https://land.copernicus.eu/pan-european/high-resolution-layers/imperviousness) (for Europe) / [National Geophysical Data Center](https://databasin.org/datasets/016d2235a5ed43ad83ceeed6c408d149) (global) and anthropogenic heat flux information from [Flanner et al. (2010)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008gl036465). All other geometrical, thermal and radiative urban canopy parameters are spatially invariant, and set to the bulk values provided in Table 1 of [Wouters et al., 2016](https://gmd.copernicus.org/articles/9/3027/2016/)).

The set of tools provided in this repo allow one to introduce LCZ-based urban canopy parameters, compiled from [Stewart and Oke (2012)](http://10.1175/BAMS-D-11-00019.1) and [Stewart et al. (2014)](http://10.1002/joc.3746).

This work is an outcome of AEVUS I and II, the COSMO Priority Tasks on "Analysis and evaluation of the TERRA_URB scheme". More info [here](http://www.cosmo-model.org/content/tasks/priorityTasks/default.htm) (project pages only accessible to COSMO members).


## Requirements
* Be a member of the [COSMO-CLM community](https://wiki.coast.hzg.de/clmcom/), in order to be able to access [EXTPAR](https://wiki.coast.hzg.de/clmcom/external-data-98599196.html).
* Have your domain file available from EXTPAR (netcdf file)
* an LCZ map covering the same region of interest. Sources for existing LCZ maps:
    * Europe: [paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0214474) | [data](http://www.wudapt.org/continental-lcz-maps/)
    * Continental United States: [paper](https://doi.org/10.1038/s41597-020-00605-z) | [data](https://doi.org/10.6084/m9.figshare.11416950.v1) 
    * An online LCZ Generator tool is currently under development; a beta version can be accessed [here](https://lcz-generator.geographie.rub.de/). Please contact Matthias.Demuzere @ rub.de for more information.
    


## Instructions

It is advised to use a python virtual environment:
1. Go into scriptdir: `cd /SCRIPT/DIR/`
2. Create virtual environment: `python3 -m venv venv` or `virtualenv venv`
3. Install module requirements: `venv/bin/pip install -r requirements.txt`
4. Use `venv/bin/pip/python` to run scripts.

The `requirements.txt` can be generated using `pipreqs`: 
```
cd /SCRIPT/DIR/
pipreqs --ignore=terra/ .
```


### Execute

The run code is currently configured for the Moscow case, as developed in Varentsov et al.
CLM and LCZ input data used in this study is provided under `data/`.

```
venv/bin/pip/python run.py
```

