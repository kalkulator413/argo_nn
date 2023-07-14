import numpy as np
import os
from netCDF4 import Dataset
from Grabber import *
from GeneralUtilities.Data.Filepath.instance import get_data_folder as get_base_folder
from datetime import datetime, timedelta, timezone

etopo = EtopoGrabber()
argo_folder = os.path.join(get_base_folder(), 'Raw/Argo')
radius = 0.5

class DataHolder:
    def __init__(self, lat:float, lon:float, temp:float, sal:float, bathymetry_grid, ssh_grid):
        self.lat = lat
        self.lon = lon
        self.temp = temp
        self.sal = sal
        self.bathymetry_grid = bathymetry_grid
        self.ssh_grid = ssh_grid

    def __str__(self):
        return f'Resurfaced at ({self.lat}, {self.lon}) with salinities {self.sal} and temperatures {self.temp}'


def loop():
    list_of_data = np.array([])

    for filename in os.listdir(argo_folder):
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

            for i in range(len(ds['POSITION_QC'])):
                # filter flag 1 for position qc
                if ds['POSITION_QC'][i] != b'1':
                    continue
                # only want real time with adjustment, not R or D
                if ds['DATA_MODE'][i] != b'A':
                    continue
                # filter flag 1 for juld qc
                if ds['JULD_QC'][i] != b'1':
                    continue

                juld = ds['JULD'][:][i]
                dt = datetime(1950, 1, 1, tzinfo=timezone.utc) + timedelta(juld)
                
                lon = ds['LONGITUDE'][:][i]
                lat = ds['LATITUDE'][:][i]

                # check if valid date and location for aviso
                try:
                    aviso = AvisoGrabber(dt.year, dt.month, dt.day)
                    aviso_grid = aviso.get_grid(lat, lon, radius)
                except:
                    continue

                # not sure how to deal with temp or psal, setting both to 0
                curr = DataHolder(lat, lon, 0, 0, aviso_grid, etopo.get_grid(lat, lon, radius))
                list_of_data = np.append(list_of_data, curr)

    return list_of_data

lst = loop()
print(f'found {len(lst)} unique data points')
