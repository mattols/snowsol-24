import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import math as m

from insolation import insolf # main package for terrain

import pygrib # NOT required
from netCDF4 import Dataset
from osgeo import gdal, osr

import random



###### GENERAL SPATIAL

def gdal_output_bounds(topo):
        geo_transform = topo.GetGeoTransform()
        return [
            geo_transform[0],
            geo_transform[3] + geo_transform[5] * topo.RasterYSize,
            geo_transform[0] + geo_transform[1] * topo.RasterXSize,
            geo_transform[3]
        ]

def hrrr_warp_domain(hrrr_file, topo_nc_path, resample_option="cubic"):
    """
    warp hrrr file based on DEM domain
    returns:
        downsampled, warped HRRR
    """
    # read in DEM data
    topo = gdal.Open(topo_nc_path, gdal.GA_ReadOnly)
    topo1 = gdal.Open(topo.GetSubDatasets()[0][0])
    spatial_info = osr.SpatialReference() # for warping
    spatial_info.ImportFromWkt(topo1.GetProjection())
    # projection info
    epsg_code = spatial_info.GetAuthorityName(None) + ":" + spatial_info.GetAuthorityCode(None)
    # temp file
    mem_hrrr_file = '/vsimem/grib_%i.tif' % random.getrandbits(32)
    # warp options for topo files
    options = gdal.WarpOptions(
                dstSRS=epsg_code,
                outputBoundsSRS=epsg_code,
                outputBounds= gdal_output_bounds(topo1),
                xRes=topo1.GetGeoTransform()[1],
                yRes=topo1.GetGeoTransform()[1],
                multithread=True,
                resampleAlg= resample_option # "bilinear"
            )
    # warp file, original, options=topo1
    gdal.Warp(mem_hrrr_file, hrrr_file, options=options)
    ## read object
    hrrr_warp = gdal.Open(mem_hrrr_file, gdal.GA_ReadOnly)
    return(hrrr_warp)
    





###### INSOLATION FUNCTIONS

def solar_positions(year, month, day, latitude=38.8697, longitude=-106.9878, 
                         timezone=-7):
    """
    solar positions fron insolf
    temporal: full day at location
    returns:
        zenith, azimuth

    default location is East River Watershed, CO
    """
    # julian day string
    jdrng = insolf.julian_day(year, month, day, np.arange(1, 25))
    # solar positions
    sunv = insolf.sunvector(jdrng, latitude, longitude, timezone)
    azimuth, zenith = insolf.sunpos(sunv)
    # normal vectors for shading (sv)
    sv = insolf.normalvector(zenith, azimuth)

    return zenith, sv, jdrng


def clearsky_insolation_inc(dem, vf, hour, jdrng, zenith, solar_vector, dlxy, visibility=60, RH=55, tempK=275.0, O3=0.02, alphag=0.2):
    """
    calculate direct and diffuse insolation on an inclined surface - clearsky model
    Temporal: single hour
    returns:
        Idir_hour, Idif_hour, I_tot (gridded)
    """
    # take hour
    current_zenith = zenith[hour]
    sv = solar_vector[:, hour]
    jday = jdrng[hour]
    
    # topographic shading effect
    sh = insolf.doshade(np.array(dem), dlxy, sv)  
    
    # gridded insolation for hour
    Idir_hour, Idiff_hour = insolf.insolation(current_zenith, jday, dem, visibility, RH, tempK, O3, alphag)
    
    # illumination angles
    hsh = insolf.hillshading(dem, dlxy, sv)
    
    # Step 6: Compute total insolation with shading (assuming Idiff_hour is not shaded)
    Idir = (Idir_hour * sh * hsh)
    Idif = (Idiff_hour * np.cos(current_zenith*m.pi/180) * vf)
    I_tot = Idir + Idif  # Diffuse insolation is not shaded

    # Return the arrays
    return Idir, Idif, I_tot







###### SPLITTING K - BOULDER

# empirical component separation

def clearness(IG, zenith, I0=1367): # Kt - clearness
    kt = IG / (I0 * m.cos(m.radians(zenith)) )
    return(kt)
               
def airmass(zenith): # m - airmass
    mass = 1/m.cos(m.radians(zenith))
    return(mass)

def localG2(kt, mass):
    a0= 0.956; a1= 1.268; a2= 3.202; a3= -6.712; a4= 2.228; a5= -0.213; a6= 0.021; 
    solve1 = (a2 + (a3*kt) + (a4*kt*kt) + (a5*mass) + (a6*mass*mass) )
    k = a0 - a1 * np.exp( -np.exp( solve1 ))
    return(k)

def diffusef_Lk2_boulder(IG, zenith_hr, I0=1367, m_elv = None):
    # clearness, airmass, and diffuse fraction (k)
    kt = clearness(IG, zenith_hr, I0)
    # airmass
    try: 
        mass = elevational_airmass(zenith_hr, m_elv)
    except:
        mass = airmass(zenith_hr)
    # diffuse fraction
    k0 = localG2(kt, mass) 
    return(k0)

def elevational_airmass(zenith, elevation):
    """
    Calculate elevation-adjusted optical air mass
    Based on Young 1987
    Uses barometric formula for pressure
    """
    # Standard air mass calculation
    base_m = 1 / m.cos(m.radians(zenith))
    
    # Elevation adjustment
    # Approximate atmospheric pressure ratio at elevation
    # Uses standard atmospheric pressure model
    # pressure_ratio = math.pow(1 - (2.25577e-5 * elevation), 5.2559)
    pressure_ratio = np.power(1 - (2.25577e-5 * elevation), 5.2559)
    
    # Adjusted air mass
    adjusted_m = base_m * pressure_ratio
    
    return adjusted_m