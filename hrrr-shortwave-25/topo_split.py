import os 
import glob
import numpy as np
from datetime import datetime, timedelta
import math as m

from netCDF4 import Dataset

from snobedo.input import HrrrParameter, SmrfTopo
from snobedo.shortwave import TopoShade
from snobedo.input import SmrfTopo

import matplotlib.pyplot as plt

from topocalc.horizon import horizon
from topocalc.shade import shade


## SPLITTING AND TOPOGRAPHIC FORCING FUNCTIONS

# FUNCTION - load terrain paramters
def static_topo_vars(ERW_TOPO_FILE):
    # static topographic variables
    with Dataset(ERW_TOPO_FILE) as dem:
        sky_view_factor = dem['sky_view_factor'][:].astype(np.float64)
        dem = dem['dem'][:].astype(np.float64)
    # topographic variables
    topo_shade = TopoShade(
        ERW_TOPO_FILE, TopoShade.SolarMethods.SKYFIELD
    )
    return dem, sky_view_factor, topo_shade

# FUNCTION - solar geom 
def hrrr_solar_geom(hrrr_dswrf, topo_shade):
    """
    Calculate solar geometry for HRRR grib file and topo shade information

    **only required for first day
    """

    hour = hrrr_dswrf.timestep_for_band(10)
    start_day = hour.strftime('%Y-%m-%d') 
    end_day = (hour + timedelta(days=1)).strftime('%Y-%m-%d')
    # hour = str(hour.strftime('%Y-%m-%d %H:%M'))
    
    # time range for solar angles
    time_range = np.arange(str(start_day), str(end_day), np.timedelta64(1, 'h'), dtype='datetime64[s]')
    
    time_range = [datetime.fromisoformat(str(r)) for r in time_range]
    topo_shade.calculate(time_range)
    
    azimuth = { key.strftime('%Y-%m-%d %H:%M'): value for key, value in topo_shade.azimuth.items() }
    zenith = { key.strftime('%Y-%m-%d %H:%M'): value for key, value in topo_shade.zenith.items() }
    incidence = { key.strftime('%Y-%m-%d %H:%M'): value for key, value in topo_shade.illumination_angles.items() }
    illumination = {
        key.strftime('%Y-%m-%d %H:%M'): value.copy() if isinstance(value, np.ndarray) else value
        for key, value in topo_shade.illumination_angles.items()
    }
    # shading corrected incidence
    illumination = horizon_shader(topo_shade, zenith, illumination) # full day - includes shading

    return azimuth, zenith, incidence, illumination

    
# FUNCTION - for shade
def horizon_shader(topo_shade, zenith, illumination_angles):
    horizon_angles = { key: horizon(value, topo_shade.topo.dem, topo_shade.topo.dx) for key, value in zenith.items() }
    # inc_copy = { key.strftime('%Y-%m-%d %H:%M'): value for key, value in illumination_angles.items() }
    for date, _val in zenith.items():
        if type(illumination_angles[date]) is np.ndarray:
            illumination_angles[date] = mask_sun(date, zenith, illumination_angles, horizon_angles) # return single date (hour)
    
    return illumination_angles

# FUNCTION - shadows from horizon - Dozier
def mask_sun(date, zenith, inc_angles, horizon_angles):
    # print(date)
    cos_z = np.cos(np.radians(zenith[date]))

    sun_angle = np.tan(np.pi / 2 - np.arccos(cos_z))
    no_sun_mask = (np.tan(np.abs(horizon_angles[date])) > sun_angle)
    inc_angles[date][no_sun_mask] = 0
    
    return inc_angles[date]

# FUNCTION - topographic splitting models (4 total)
def toposplitModels4(hrrr_dswrf, zenith, azimuth, incidence_angles, illumination, sky_view_factor, topo_shade):
    """
    Topo models for HRRR
        Provides splitting method based on 

    Returns:
        ghi  - global horizontal radiation (flat model
        dsw1 - + illumination angle (inc)
        dsw2 - + inc + skyview (vf)
        dsw3 - + inc + vf + shading (COMPLETE)
        
    """
    # define hour
    hour = hrrr_dswrf.timestep_for_band(10)
    hour = str(hour.strftime('%Y-%m-%d %H:%M'))
    
    # variables
    ghi = hrrr_dswrf.grib_file.GetRasterBand(10).ReadAsArray()
    k = hrrr_dswrf.grib_file.GetRasterBand(11).ReadAsArray() / hrrr_dswrf.grib_file.GetRasterBand(12).ReadAsArray()
    
    ## Topo models
    # DSW1 - No cast shadows - no view factor
    dni0 = ((ghi * (1-k)) / (np.cos(zenith[hour]*m.pi/180))) * incidence_angles[hour]
    dif0 = (dni0 * k) 
    dsw1 = dni0 + dif0
    # DSW2 - No cast shadows
    dif = dif0 * sky_view_factor
    dsw2 = dni0 + dif
    # DSW3 - Final model
    dni = ((ghi * (1-k)) / (np.cos(zenith[hour]*m.pi/180))) * illumination[hour]
    dsw3 = dni + dif

    return ghi, dsw1, dsw2, dsw3
    

# # FUNCTION - Wrapper for single day
# def toposplit_day(hrrr_dswrf):
#     """
#     Wrapper for toposplit over single day _ NEEDS TO ACCOUNT FOR UTC CROSSOVER TO NEW DAY
    
#     """
#     # check topo vars
#     if 'topo_shade' not in globals():      # if 'my_variable' not in locals():
#         dem, sky_view_factor, topo_shade = static_topo_vars(ERW_TOPO_FILE)

#     # calculate solar geometry (FUll day)
#     azimuth, zenith, incidence = hrrr_solar_geom(hrrr_dswrf, topo_shade)

#     # run for all days
#     for i in range(zenith):
#         ghi, dsw1, dsw2, dsw3 = toposplitModels4(hrrr_dswrf, zenith, azimuth, illumination, sky_view_factor, topo_shade)



# FUNCTION - wrapper for single day - NEEDS TO ACCOUNT FOR PREVIOUS DAY DAYLIGHT HOURS AT early UTC time
def toposplit_day_wrap_to_netcdf(ERW_TOPO_FILE, HRRR_DIR, output_root_dir, resample_method='cubic'):
    """
    Runs toposplit on HRRR files for a full day (UTC) and save each of ghi, dsw1, dsw2, dsw3 as separate NetCDF files
    under different subfolders.

    Args:
        ERW_TOPO_FILE (str): Path to topo.nc
        HRRR_DIR (str): Path to HRRR folder (e.g., .../hrrr.20220301/)
        output_root_dir (str): Root directory where outputs will be stored
    """
    # bug fix
    os.environ['PROJ_LIB'] = '/uufs/chpc.utah.edu/common/home/skiles-group1/jmeyer/conda/envs/smeshr/share/proj'

    # extract date
    date_str = os.path.basename(HRRR_DIR).split('.')[-1]  # '20220301'
    base_date = datetime.strptime(date_str, "%Y%m%d")
    
    # static terrain data
    if 'topo_shade' not in globals() or 'topo_shade' not in locals():
        dem, sky_view_factor, topo_shade = static_topo_vars(ERW_TOPO_FILE)

    # list of HRRR files for the day (UTC forecast 6)
    hrrr_files = sorted(glob.glob(os.path.join(HRRR_DIR, 'hrrr.t*z.wrfsfcf06.grib2')))
    if not hrrr_files:
        raise FileNotFoundError(f"No HRRR files found in {HRRR_DIR}")

    # solar geometry once for day
    sample_dswrf = HrrrParameter(ERW_TOPO_FILE, hrrr_files[0], 'cubic')
    azimuth, zenith, incidence_angles, illumination = hrrr_solar_geom(sample_dswrf, topo_shade)

    # shape
    sample_band = sample_dswrf.grib_file.GetRasterBand(10).ReadAsArray()
    n_y, n_x = sample_band.shape
    n_t = len(hrrr_files)

    # time series arrays
    ghi_stack = np.zeros((n_t, n_y, n_x), dtype=np.float32)
    dsw1_stack = np.zeros_like(ghi_stack)
    dsw2_stack = np.zeros_like(ghi_stack)
    dsw3_stack = np.zeros_like(ghi_stack)
    time_list = []

    for i, file_path in enumerate(hrrr_files):
        try:
            hrrr_dswrf = HrrrParameter(ERW_TOPO_FILE, file_path, resample_method)
            hour = hrrr_dswrf.timestep_for_band(10)
            time_list.append(hour)

            ghi, dsw1, dsw2, dsw3 = toposplitModels4(
                hrrr_dswrf, zenith, azimuth, incidence_angles, illumination, sky_view_factor, topo_shade
            )

            ghi_stack[i] = ghi
            dsw1_stack[i] = dsw1
            dsw2_stack[i] = dsw2
            dsw3_stack[i] = dsw3

            print(f"[✓] Processed {hour.strftime('%Y-%m-%d %H:%M')}")

        except Exception as e:
            print(f"[!] Error processing {file_path}: {e}")

    # create and write each NetCDF file
    def save_stack_to_netcdf(var_name, data_stack, out_dir):
        os.makedirs(out_dir, exist_ok=True)
        output_file = os.path.join(out_dir, f'{var_name}_{date_str}.nc')

        with Dataset(output_file, 'w', format='NETCDF4') as nc:
            # dimensions
            nc.createDimension('time', n_t)
            nc.createDimension('y', n_y)
            nc.createDimension('x', n_x)

            # time
            time_var = nc.createVariable('time', 'f8', ('time',))
            time_var.units = f'hours since {date_str[:4]}-{date_str[4:6]}-{date_str[6:]} 00:00:00'
            time_var.calendar = 'standard'
            time_var[:] = [(t - time_list[0]).total_seconds() / 3600.0 for t in time_list]

            # data variable
            var = nc.createVariable(var_name, 'f4', ('time', 'y', 'x'), zlib=True, complevel=4)
            var.units = 'W/m^2'
            var[:, :, :] = data_stack

            # metadata
            nc.description = f'{var_name} output from HRRR topographic split model'
            nc.history = f'Created {datetime.now()}'
            nc.source = 'HRRR model + custom topo-splitting'

        print(f"[✓] Saved: {output_file}")

    # save each stack in its own folder
    save_stack_to_netcdf('ghi', ghi_stack, os.path.join(output_root_dir, 'hrrr_dswrf'))
    save_stack_to_netcdf('dsw1', dsw1_stack, os.path.join(output_root_dir, 'hrrr_dsw1'))
    save_stack_to_netcdf('dsw2', dsw2_stack, os.path.join(output_root_dir, 'hrrr_dsw2'))
    save_stack_to_netcdf('dsw3', dsw3_stack, os.path.join(output_root_dir, 'hrrr_dsw3'))

    print("\n All NetCDF files saved successfully.")

