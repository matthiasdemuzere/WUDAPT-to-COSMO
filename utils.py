import pandas as pd
import numpy as np
import xarray as xr
from scipy.interpolate import RegularGridInterpolator
import scipy.ndimage

## Helper function to prepare the urban canopy data.
def prepareUCPlookup(ucpFile, saiWeight=False, snow_f=0, alb_snow=0.70, emi_snow=0.997):

    """

        AUTHOR: Matthias Demuzere (matthias.demuzere [@] rub.de)

        VERSION: 1.0

        DATE: 2020-07-01

        :param ucpFile  : Absolute path to urban canopy parameter csv file.
        :param saiWeight: Weigh parameters according to Surface Area Index (Default = False)
        :param snow_f   : snow fraction (default = 0)
        :param alb_snow : snow albedo (default = 0.7)
        :param emi_snow : emissivity albedo (default = 0.997)
        :return:

        :INFO:
        Code generally follows SURY: https://github.com/hendrikwout/sury/blob/master/sury.py

        * Wouters, H., Demuzere, M., Blahak, U., Fortuniak, K., Maiheu., B.,
            Camps, J., Tielemans, and N. P. M. van Lipzig, 2016.  The efficient
            urban-canopy dependency parametrization SURY (v1.0) for atmospheric modelling:
            description and application with the COSMO-CLM model (v5.0_clm6) for a
            Belgian Summer, Geosci. Model Dev., 2016.

        Define the look-up table, based on the values of:
        * Stewart, I. D., & Oke, T. R. (2012). Local Climate Zones for Urban Temperature Studies.
            Bulletin of the American Meteorological Society, 93(12), 1879–1900.
        * Stewart, I. D., Oke, T. R., & Krayenhoff, E. S. (2014). Evaluation of the ‘local climate zone’
            scheme using temperature observations and model simulations. International Journal of Climatology,
            34(4), 1062–1080. https://doi.org/10.1002/joc.3746

        The latter paper describes thermal admittance values of facets only. Heat conductivity and
        capacity values are obtained via Scott Krayenhoff (personal communication).

    """

    ## Read look-up table
    ucp = pd.read_csv(ucpFile,sep=';',index_col=0).iloc[:17,:]

    ## Canyon albedo reduction factor, eq. 15 Wouters et al. (2016)
    psi_canyon = np.exp(-0.6 * ucp['URB_H2W'])
    psi_canyon[10:] = 0 # Set to zero for non-urban LCZs

    ## Total albedo reduction factor, eq. 14 Wouters et al. (2016)
    #psi_bulk = ucp['URB_BLDFR'] + (1-ucp['URB_BLDFR'])*psi_canyon*ucp['URB_H2W']

    ## Bulk albedo
    alb_roof_snow = ucp['URB_RfALB'] * (1. - snow_f) + alb_snow * snow_f
    alb_road_snow = ucp['URB_RdALB'] * (1. - snow_f) + alb_snow * snow_f
    alb_wall_snow = ucp['URB_WaALB'] * (1. - snow_f) + alb_snow * snow_f
    ucp['URB_SALB'] = (alb_road_snow + 2. * ucp['URB_H2W'] * alb_wall_snow) / \
                            (1. + 2. * ucp['URB_H2W']) * psi_canyon * (1. - ucp['URB_BLDFR']) \
                            + alb_roof_snow * ucp['URB_BLDFR']
    ucp.loc[11:,'URB_SALB'] = 0
    ucp['URB_TALB'] = ucp['URB_SALB'].copy()

    ## Bulk emissivity
    emi_roof_snow = (1. - ucp['URB_RfEMI']) * (1. - snow_f) + (1. - emi_snow) * snow_f
    emi_road_snow = (1. - ucp['URB_RdEMI']) * (1. - snow_f) + (1. - emi_snow) * snow_f
    emi_wall_snow = (1. - ucp['URB_WaEMI']) * (1. - snow_f) + (1. - emi_snow) * snow_f
    ucp['URB_EMIS'] = 1. - ((emi_road_snow + 2. * ucp['URB_H2W'] * emi_wall_snow) \
                     / (1. + 2. * ucp['URB_H2W']) * psi_canyon * (1. - ucp['URB_BLDFR']) \
                     + emi_roof_snow * ucp['URB_BLDFR'])
    ucp.loc[11:,'URB_EMIS'] = 0


    ## Calculate Surface Area Index from geometrical considerations (Eq. 3)
    SAI = (1. + 2. * ucp['URB_H2W']) * (1. - ucp['URB_BLDFR']) + ucp['URB_BLDFR']

    ## Get bulk Heat capacity and conductivity, using eq. 10, 11 and 4.
    ucp['URB_HCON'] = ((1-ucp['URB_BLDFR']) / SAI ) * \
              (2*ucp['URB_H2W']*ucp['URB_WaHCON'] + ucp['URB_RdHCON']) + \
              ( ucp['URB_BLDFR'] / SAI * ucp['URB_RfHCON'])
    ucp['URB_HCAP'] = ((1-ucp['URB_BLDFR']) / SAI ) * \
              (2*ucp['URB_H2W']*ucp['URB_WaHCAP'] + ucp['URB_RdHCAP']) + \
              ( ucp['URB_BLDFR'] / SAI * ucp['URB_RfHCAP'])

    ## Also add the thermal admittance
    #ucp['URB_TADM'] = (ucp['URB_HCAP']*ucp['URB_HCON'])**0.5

    ## iS SAI weighting requested, according to Eq. 4?
    ## This is done within TERRA_URB, so no need to do for COSMO/CLM input files.
    if saiWeight:
        ucp['URB_HCON'] = ucp['URB_HCON'] * SAI
        ucp['URB_HCAP'] = ucp['URB_HCAP'] * SAI
        #ucp['URB_TADM'] = ucp['URB_TADM'] * SAI

    return ucp


## Helper function to do the interpolation
def cosmoInterpolator(xLcz, yLcz, dataLcz, xClm, yClm, interpMethod='linear',
                      aggregation=True, aggregationScale = 2):

    """
        AUTHOR: Matthias Demuzere (matthias.demuzere [@] rub.de)

        VERSION: 1.0

        DATE: 2020-07-01

        :param xLcz: values of LCZ longitudes
        :param yLcz: values of LCZ latitudes
        :param dataLcz: 2D array of LCZ parameter value
        :param xClm: values of COSMO-CLM longitudes
        :param yClm: values of COSMO-CLM latitudes
        :param interpMethod: "linear" (default) or "nearest"
        :param aggregation: Boolean to aggregate or not
        :param aggregationScale: scaler integer to do a pre-processing aggregation (default = 2, ~ half size of CLM dimensions)

        :return: 2d-array of the lcz paramater on the COSMO-CLM grid / dimensions
    """

    if aggregation:
        ## First aggregate the LCZ data, to ~ twice the size of the CLM domain (defined by aggregationParam)
        aggNr = np.round(np.max([len(xClm) / len(xLcz), len(yClm) / len(yLcz)]), 3)
        lczVarAgg = scipy.ndimage.zoom(dataLcz, aggNr * aggregationScale, order=1);

        ## Get corresponding new LCZ lat and lon values
        latsAgg = np.linspace(yLcz.min(), yLcz.max(), lczVarAgg.shape[0])
        lonsAgg = np.linspace(xLcz.min(), xLcz.max(), lczVarAgg.shape[1])

        interp_object = RegularGridInterpolator((latsAgg, lonsAgg), lczVarAgg, method=interpMethod)

    else:
        interp_object = RegularGridInterpolator((yLcz, xLcz), dataLcz, method=interpMethod)

    return interp_object, aggNr



def lcz_to_cosmo(ucpFile, clmFile, lczFile, bandNr, ucpVersion, nrLcz=17,
                 interpMethod='linear', aggregation=True, aggregationScale=2,
                 isaWeight=True, saiWeight=False,
                 fileNameExt=''):
    """
        Function to introduce LCZ Urban Canopy Parameters into CCLM domain file

        AUTHOR: Matthias Demuzere (matthias.demuzere [@] rub.de)

        VERSION: 1.0

        DATE: 2020-07-01

        :param ucpFile: full absolute path name to ucp .csv table.
        :param clmFile: full absolute path name to COSMO-CLM domain file
        :param lczFile: full absolute path name to lcz geotiff file.
        :param bandNr: integer, referring to version of LCZ map:
                0 = lcz, 1 = lczFilter, 2 = lczFilter CGLS mask, 3 = lczFilter GAIA mask
        :param ucpVersion: version of ucp file used: 'high' (Stewart and Oke, 2014) or 'default'
        :param gcFile: full path to file containing Globcover parameters per class
        :param nrLcz: highest value of LCZ class present (default is 17)
        :param interpMethod: 'linear' (default) or 'nearest'
        :param aggregation: Boolean to aggregate or not
        :param aggregationScale: scaler integer to do a pre-processing aggregation (default = 2, ~ half size of CLM dimensions)
        :param isaWeight: Boolean. Weighs parameter according to ISA fraction (default == True)
        :param saiWeight: Weigh parameters according to Surface Area Index (Default = False)
        :param fileNameExt: provide opportunity for additional file name extension (default: '')

        :return:
        LAFxxxxxx_lcz.nc file, with
            - additional fields for 'ISA','AHF','BLDH','BLDFR','HW','CVS','ALB'
            - correcting for double counting for Globcover affected fields

        References:
        - Stewart, I.D., Oke, T.R., 2012. Local Climate Zones for Urban Temperature Studies.
            Bull. Am. Meteorol. Soc. 93, 1879–1900. https://doi.org/10.1175/BAMS-D-11-00019.1
        - Wouters, H., Demuzere, M., Blahak, U., Fortuniak, K., Maiheu, B., Camps, J.,
            Tielemans, D., van Lipzig, N.P.M., 2016. Efficient urban canopy parametrization
            for atmospheric modelling: description and application with the COSMO-CLM model
            (version 5.0_clm6) for a Belgian Summer. Geosci. Model Dev. 9, 3027–3054.
            https://doi.org/10.5194/gmd-2016-58

    """

    ## for testing
    #nrLcz = 17; interpMethod = 'linear'; aggregation = True; aggregationScale = 2; isaWeight = True

    lookupUCP  = prepareUCPlookup(ucpFile,saiWeight)

    ## Read lcz file, make copy of original domainFile
    lczMap = xr.open_rasterio(lczFile)[bandNr,:,:].astype('int')
    lczMap = lczMap.rename({'x': 'lon', 'y': 'lat'})
    lczMap = lczMap.reindex(lat=lczMap.lat[::-1])

    ## Read LCZ map coordinates
    xLcz, yLcz = lczMap.lon.values, lczMap.lat.values

    ## Create a new domain file as a copy of the original one
    clmFileNew = clmFile.replace('.nc','_lcz_{}_{}{}.nc'.format(bandNr,ucpVersion,fileNameExt))
    clm  = xr.open_dataset(clmFile)

    ## Read COSMO-CLM coordinates
    xClm, yClm = clm.lon.values, clm.lat.values

    ## Define list of all urban parameters.
    ## Only change FR_PAVED and URBAN after fixing double counting.
    urbParameters = ['ISA',
                     'FR_PAVED',
                     'URB_BLDFR',
                     'URB_BLDH',
                     'URB_H2W',
                     'AHF',
                     'URB_SALB',
                     'URB_TALB',
                     'URB_EMIS',
                     'URB_HCON',
                     'URB_HCAP']

    ## Create maps of UCPs, store in xarray data object
    for ucp in urbParameters:

        ## If variables not present, create variable with dummy values
        ## in dataarray, and overwrite below
        if not ucp in list(clm.variables.keys()):
            clm[ucp] = clm['URBAN'].copy()

        keys = np.arange(1, nrLcz+1, 1)
        values = lookupUCP[ucp]
        out = np.empty((max(keys) + 1,), object); out[list(keys)] = values
        dataLcz = np.array(out[lczMap], dtype='float')

        ## Set nans to 0, required for interpolation.
        dataLcz[np.isnan(dataLcz)] = 0

        ## Store isa separately for weighting
        if ucp == 'ISA':
            dataISA = dataLcz

        ## Weigh values according to ISA -  Step 1
        if isaWeight == True and not ucp in ['ISA','AHF','FR_PAVED']:
            print("isaWeight on: ucp's are being weighed by ISA fraction")
            dataLcz = dataLcz * dataISA

        ## Get interpolation object
        interp_object, aggNr  = cosmoInterpolator(xLcz, yLcz, dataLcz, xClm, yClm,
                                           interpMethod, aggregation, aggregationScale)

        ## Apply to get resampled data
        clmPoints = yClm, xClm
        dataLczResampled = interp_object(clmPoints)

        ## Store isa lczClm separately for weighting
        if ucp == 'ISA':
            dataISAres = dataLczResampled

        ## re-Weigh values according to ISA -  Step 2
        if isaWeight == True and not ucp in ['ISA','AHF','FR_PAVED']:
            dataLczResampled = dataLczResampled / dataISAres

        ## Set all nans to zero, otherwise issues with COSMO.
        dataLczResampled[np.isnan(dataLczResampled)] = 0

        ## Add values to file.
        clm[ucp].values = dataLczResampled

        print('{} field has been updated'.format(ucp))

        ## Change variable attributes
        clm[ucp].attrs['data_set'] = 'Values derived from Local Climate Zone properties'

    ## Set FR_PAVED equal to ISA.
    clm['FR_PAVED'].values = clm['ISA'].values

    ## Add global attribute
    clm.attrs['note2'] = 'LCZ Urban canopy look up data retrieved from Stewart and Oke (2012) \n' \
                            ' and Stewart et al. (2014). Conversion to bulk properties done via SURY \n' \
                            'from Wouters et al. (2016): https://github.com/hendrikwout/sury/blob/master/sury.py.'

    ## Write to file
    print("Writing COSMO-CLM domain with LCZ values to {}".format(clmFileNew))
    clm.to_netcdf('{}'.format(clmFileNew))

    return '{}'.format(clmFileNew)


def removeDoubleCounting(clmFile,gcFile,removeUrban=True,qLow=0.25,qHigh=0.75,fileNameExt=''):

    """
        Function to remove the double counting of URBAN-BASED parameter values

        AUTHOR: Matthias Demuzere (matthias.demuzere [@] rub.de)

        VERSION: 1.0

        DATE: 2020-07-01

        :param clmFile: full absolute path name to COSMO-CLM domain file
        :param gcFile: full absolute path name to globcover look-up table
        :param removeUrban: Boolean to indicate if URBAN effect from GlobCover needs to be removed (default=True)
                            If False, the procedure from EXTPAR is reconstructed, values should =~ input file.
        :param qLow, qHigh: Low and high quantile.
                            Where URBAN == 1, values are random sampled between qLow qnd qHigh (from URBAN == 0 pixels)
        :param fileNameExt: provide opportunity for additional file name extension (default: '')

        :return:
        Adjusted COSMO/CLM domain file, with:
        * fixes for (if enabeld): 'Z0', 'PLCOV_MN', 'PLCOV_MX', 'LAI_MN', 'LAI_MX', 'ROOTDP',
                                  'EMIS_RAD', 'SKC', 'RSMIN', 'FOR_D', 'FOR_E'
        * URBAN and FR_PAVED set to ISA.
    """

    ## Constants
    hp = 30 # height of Prandtl-layer, taken from EXTPAR's mo_agg_globcover.f90

    ## Read original file, adjusted for LCZs
    clm  = xr.open_dataset('{}'.format(clmFile)) #, decode_coords=False)

    ## Make string for output file, allow for file name extensions
    clm_oFile = clmFile.replace('.nc', '_fixDC_{}{}.nc'.format(removeUrban,fileNameExt))

    ## Read relevant fields for double counting
    lu = clm.LU_CLASS_FRACTION.values
    urb = clm.URBAN.values
    frl = clm.FR_LAND.values

    ## Read parameterfile
    gc_lookup = pd.read_csv(gcFile, sep=';').iloc[:,2:]
    gc_vars = gc_lookup.columns.to_list()

    ## Whether or not to fix the double counting
    if removeUrban:
        print('Double counting will be removed from the urban pixels')
        nonUrban = [x for x in range(23) if x != 18]
    else:
        print('Double counting not addressed: reconstruction of original domain file values.')
        nonUrban = list(range(23))

    ## Define the pixels that need to be altered
    touchPix = np.asarray(np.logical_and(np.logical_and(urb > 0, urb != 1), frl > 0.5))
    print('{} pixels identified with 0 < urb < 1, frl > 0.5'.format(np.sum(touchPix)))

    ## Start replacing the values, use array broadcasting for efficiency
    for gc_var in gc_vars:
        if gc_var in list(clm.var()):
            print('Fixing double counting for {}: {}'.format(gc_var, removeUrban))

            ## Initialize array to store new values in.
            clmValue = clm[gc_var].values
            tmp = clmValue.copy()

            #tmp[np.isnan(clmValue)] = np.nan

            ## Stretch land use to unity if urban fractions are removed.
            luStretched = lu[nonUrban, :, :] * (1 / np.sum(lu[nonUrban, :, :], axis=0))

            if gc_var in ['PLCOV_MX', 'PLCOV_MN', 'EMIS_RAD', 'SKC']:
                tmp_v = np.sum(luStretched * np.expand_dims(gc_lookup[gc_var][nonUrban], axis=[1, 2]), axis=0)
                tmp[touchPix] = tmp_v[touchPix]

            elif gc_var in ['Z0']:
                tmp_v = np.sum(
                    luStretched / np.expand_dims(np.log(hp) - np.log(gc_lookup[gc_var][nonUrban]), axis=[1, 2]),
                    axis=0)
                tmp[touchPix] = hp * np.exp(-1 / tmp_v[touchPix])

            elif gc_var in ['ROOTDP', 'LAI_MN', 'LAI_MX', 'FOR_D', 'FOR_E']:
                luStretched = lu[nonUrban, :, :] * (1 / np.sum(lu[nonUrban, :, :]))
                tmp_n = np.sum(
                    luStretched * np.expand_dims(gc_lookup['PLCOV_MX'][nonUrban] * gc_lookup[gc_var][nonUrban], axis=[1, 2]),
                    axis=0)
                tmp_d = np.sum(luStretched * np.expand_dims(gc_lookup['PLCOV_MX'][nonUrban], axis=[1, 2]), axis=0)
                tmp[touchPix] = tmp_n[touchPix] / tmp_d[touchPix]

            ## Fix the pixels with URB == 1, replace random values between Q1 and Q2 over domain (non-urban)
            if removeUrban:
                replacePixels = np.logical_and(urb == 0, lu[-1, :, :] == 0)
                #tmp[urb == 1] = np.nanmedian(tmp[replacePixels])
                tmp_q = np.quantile(tmp[replacePixels], [qLow,qHigh])
                tmp[urb == 1] = np.random.uniform(tmp_q[0],tmp_q[1],size=tmp.shape)[urb == 1]

            ## Set values in clm file
            clm[gc_var].values = tmp
            print('Done for {}'.format(gc_var))

    ## Set URBAN and FR_PAVED to ISA, for consistency
    clm['URBAN'].values = clm['ISA'].values
    clm['FR_PAVED'] = clm['ISA'].copy()

    ## Write to file
    print("Writing COSMO-CLM domain with double counting fixed: {}".format(clm_oFile))
    clm.to_netcdf('{}'.format(clm_oFile))
