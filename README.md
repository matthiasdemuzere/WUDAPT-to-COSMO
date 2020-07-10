# LCZ-to-CCLM

Set of tools to use LCZ-based urban canopy parameters in DWD's COSMO-CLM NWP and regional climate model.

## Context
TERRA_URB, the urban canopy parameterization in COSMO-CLM ([Wouters et al., 2016](https://gmd.copernicus.org/articles/9/3027/2016/)), by default uses impervious surface area information from ESA (Europe) / NOAA (global) and anthropogenic heat flux information from [Flanner et al. (2010)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008gl036465). All other geometrical, thermal and radiative urban canopy parameters are spatially invariant, and set to the values provided in Table of 1 [Wouters et al., 2016](https://gmd.copernicus.org/articles/9/3027/2016/)).

