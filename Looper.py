import numpy as np
import os
from netCDF4 import Dataset
from argo_nn.Grabber import *
from argo_nn.TS_Parser import *
from GeneralUtilities.Data.Filepath.instance import get_data_folder as get_base_folder
from datetime import datetime, timedelta, timezone
import pickle

etopo = EtopoGrabber()
argo_folder = os.path.join(get_base_folder(), 'Raw/Argo')
radius = 3
pres = sal_ds['PRESSURE'][pressure_range]

class DataHolder:
    def __init__(self, temp_sal, bathymetry_grid, ssh_grid, angle:float, magnitude:float):
        self.temp_sal = temp_sal
        self.bathymetry_grid = bathymetry_grid
        self.ssh_grid = ssh_grid
        self.angle = angle
        self.magnitude = magnitude

def preloop():
    # read in EOFs for temperature and salinity data
    # and store the first 10 of each
    with open('eofs/u.pkl', 'rb') as f:
        eofs = pickle.load(f).T

    return eofs[:6]
    # add any further optimizations ...

def interpolate(argo_depths, rg_depths, argo_readings, pres_qc):
    result = np.zeros(len(rg_depths))
    prev = 0
    for i, depth in enumerate(rg_depths):

        while prev < len(argo_depths) - 1 and not (argo_depths[prev] < depth and argo_depths[prev+1] >= depth):
            if (not type(argo_depths.mask) == np.bool_) and argo_depths.mask[prev+1]:
                print(f'couldnt exceed depth {depth}')
                return None
            prev += 1

        if (pres_qc[prev] != b'1') or (pres_qc[prev+1] != b'1'):
                print(f'failed quality control')
                return None
        x1 = argo_depths[prev]
        x2 = argo_depths[prev+1]
        y1 = argo_readings[prev]
        y2 = argo_readings[prev+1]

        slope = (y2 - y1) / (x2 - x1)
        interpolated_val = slope * (depth - x1) + y1
        result[i] = interpolated_val

    return result

def loop():
    list_of_data = np.array([])
    eofs = preloop()

    for filename in os.listdir(argo_folder)[:1]:
        subfolder = os.path.join(argo_folder, filename)
        print('Currently traversing', filename)

        # for the sake of testing, only go through 1 file in each subfolder
        for f in os.listdir(subfolder)[:1]:
            print('Looking at float', f)

            prof_file = os.path.join(subfolder, f, f'{f}_prof.nc')
            if not os.path.isfile(prof_file):
                # print(prof_file)
                continue
            ds = Dataset(prof_file)

            for i in range(len(ds['POSITION_QC']) - 1):
                # filter flag 1 for position qc
                if ds['POSITION_QC'][i] != b'1':
                    continue
                # reject real time data
                if ds['DATA_MODE'][i] == b'R':
                    continue
                # filter flag 1 for juld qc
                if ds['JULD_QC'][i] != b'1':
                    continue

                # check if next resurface height and time good
                if ds['POSITION_QC'][i+1] != b'1':
                    continue
                if ds['JULD_QC'][i+1] != b'1':
                    continue

                juld = ds['JULD'][:][i]
                juld_next = ds['JULD'][:][i+1]
                if (juld_next - juld) < 9 or (juld_next - juld) > 11:
                    continue
                
                float_pres = ds['PRES'][i]
                temp = interpolate(float_pres, pres, ds['TEMP'][i], ds['PRES_QC'][i])
                sal = interpolate(float_pres, pres, ds['PSAL'][i], ds['PRES_QC'][i])
                if sal is None or temp is None:
                    continue

                dt = datetime(1950, 1, 1, tzinfo=timezone.utc) + timedelta(juld)
                
                lon = ds['LONGITUDE'][:][i]
                lat = ds['LATITUDE'][:][i]

                try:
                    aviso_grid = AvisoGrabber(dt.year, dt.month, dt.day).get_grid(lat, lon, radius)
                    etopo_grid = etopo.get_grid(lat, lon, radius)
                except Exception as e:
                    # print(e)
                    continue

                dlon = ds['LONGITUDE'][:][i+1] - lon
                dlat = ds['LATITUDE'][:][i+1] - lat

                # longitude rollover
                if dlon > 100:
                    continue

                angle = np.arctan(dlat / dlon)

                # quadrant 2 fix
                if (dlat >= 0 and dlon < 0):
                    angle += np.pi
                # quadrant 3 fix
                if (dlat < 0 and dlon < 0):
                    angle -= np.pi

                magnitude = np.sqrt(dlat**2 + dlon**2)

                ts = np.append(temp, sal)
                proj_ts = [np.dot(x, ts) for x in eofs]

                curr = DataHolder(proj_ts, aviso_grid, etopo_grid, angle, magnitude)
                list_of_data = np.append(list_of_data, curr)

    return list_of_data

lst = loop()
print(f'found {len(lst)} unique data points')
