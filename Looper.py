import numpy as np
import os
from netCDF4 import Dataset
from argo_nn.Grabber import *
from argo_nn.TS_Parser import *
from GeneralUtilities.Data.Filepath.instance import get_data_folder as get_base_folder
from datetime import datetime, timedelta, timezone

etopo = EtopoGrabber()
argo_folder = os.path.join(get_base_folder(), 'Raw/Argo')
radius = 0.5

# TODO:
# make etopo grabber faster by having a lookup table for all possible (lat, lon)
# make aviso grabber faster by reading in all possible (year, month) combinations beforehand

class DataHolder:
    def __init__(self, lat:float, lon:float, temp, sal, dlat:float, dlon:float, bathymetry_grid, ssh_grid, angle:float, magnitude:float):
        self.lat = lat
        self.lon = lon
        self.temp = temp
        self.sal = sal
        self.dlat = dlat
        self.dlon = dlon
        self.bathymetry_grid = bathymetry_grid
        self.ssh_grid = ssh_grid

    def __str__(self):
        return f'Resurfaced at ({self.lat}, {self.lon}) with salinities {self.sal} and temperatures {self.temp}'


def loop():
    list_of_data = np.array([])

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
                # no real time data
                if ds['DATA_MODE'][i] == b'R':
                    continue
                # filter flag 1 for juld qc
                if ds['JULD_QC'][i] != b'1':
                    continue

                # check if next resurface height good
                if ds['POSITION_QC'][i+1] != b'1':
                    continue

                # check if valid date and location for aviso

                juld = ds['JULD'][:][i]
                juld_next = ds['JULD'][:][i+1]
                if (juld_next - juld) < 9 or (juld_next - juld) > 11:
                    continue

                dt = datetime(1950, 1, 1, tzinfo=timezone.utc) + timedelta(juld)
                
                lon = ds['LONGITUDE'][:][i]
                lat = ds['LATITUDE'][:][i]

                try:
                    aviso_grid = AvisoGrabber(dt.year, dt.month, dt.day).get_grid(lat, lon, radius)
                    temp, sal = TS_Grabber(dt.year, dt.month).get_profiles(lat, lon)
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

                # not sure how to deal with temp or psal, setting both to 0
                curr = DataHolder(lat, lon, temp, sal, dlat, dlon, aviso_grid, etopo.get_grid(lat, lon, radius), angle, magnitude)
                list_of_data = np.append(list_of_data, curr)

    return list_of_data

lst = loop()
print(f'found {len(lst)} unique data points')
