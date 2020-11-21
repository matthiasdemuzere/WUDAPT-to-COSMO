import os
abspath = os.path.abspath(__file__)
WORKDIR = os.path.dirname(abspath)
os.chdir(WORKDIR)
from utils import lcz_to_cosmo, remove_double_counting

## Adjust directory here
CLMFILE = f"{WORKDIR}/data/MSK_0.009bg3_Globcover.nc"
LCZFILE = f"{WORKDIR}/data/LCZ_Russia_Moscow.tif"

# FILES
UCPFILE = f"{WORKDIR}/tables/LCZ_UCP_default.csv"
GCFILE  = f"{WORKDIR}/tables/globcover_lookup.csv"


# EXECUTE FUNCTIONS
# 1. Assign UCP values to LCZ map and convert to COSMO Grid
CLM_FILE_NEW = lcz_to_cosmo(
    ucpFile=UCPFILE,
    clmFile=CLMFILE,
    lczFile=LCZFILE,
    bandNr=3,
    ucpVersion='default',
    nrLcz=17,
    interpMethod='linear',
    aggregation=True,
    aggregationScale=2,
    isaWeight=True,
    saiWeight=False,
    fileNameExt='_Varentsov_etal_Atm')

# 2. Address the double counting issue.
remove_double_counting(
    clmFile=CLM_FILE_NEW,
    gcFile=GCFILE,
    removeUrban=True)