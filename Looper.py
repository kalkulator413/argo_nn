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
num_eofs = 8

def preloop():
    with open('eofs/u.pkl', 'rb') as f:
        eofs = pickle.load(f).T

    with open('eofs/means.pkl', 'rb') as f:
        means = pickle.load(f)

    return eofs[:num_eofs], means
    # add any further optimizations ...

def loop():
    list_of_data = np.array([])
    eofs, means = preloop()

    for fnum, filename in enumerate(os.listdir(argo_folder)):
        subfolder = os.path.join(argo_folder, filename)
        print('Currently traversing', filename)

        # for the sake of testing, only go through 1 file in each subfolder
        for f in os.listdir(subfolder)[:1]:
        # for f in os.listdir(subfolder):
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
                float_pres_qc = ds['PRES_QC'][i]

                t = interpolate(float_pres, pres, ds['TEMP'][i], float_pres_qc, ds['TEMP_QC'][i])
                s = interpolate(float_pres, pres, ds['PSAL'][i], float_pres_qc, ds['PSAL_QC'][i])

                if t is None or s is None:
                    continue
                
                concat = np.append(t / np.sqrt(tempvar), s / np.sqrt(salvar)) - means
                ts = np.array([np.dot(concat, eofs[i]) for i in range(num_eofs)])

                dt = datetime(1950, 1, 1, tzinfo=timezone.utc) + timedelta(juld)
                
                lon = ds['LONGITUDE'][:][i]
                lat = ds['LATITUDE'][:][i]

                try:
                    # aviso_grid = AvisoGrabber(dt.year, dt.month, dt.day).get_grid(lat, lon, radius)
                    # etopo_grid = etopo.get_grid(lat, lon, radius)
                    pass
                except Exception as e:
                    # print(e)
                    continue

                # longitude rollover
                if abs(ds['LONGITUDE'][:][i+1] - lon) > 100:
                    continue

                ssh_slope = 0
                ssh_dir = 0
                bath_slope = 0
                bath_dir = 0
                roughness = 0
                curr = np.append(ts, np.array([ssh_slope, ssh_dir, bath_slope, bath_dir, roughness]))
                curr = np.append([fnum,int(f),i,dt.year,dt.month,dt.day,lat,lon,ds['LATITUDE'][:][i+1], ds['LONGITUDE'][:][i+1]], curr)

                list_of_data = np.append(list_of_data, curr)

    list_of_data = np.reshape(list_of_data, (-1, num_eofs + 5 + 10))
    list_of_data[:,[0,1]] = list_of_data[:,[0,1]].astype(int)
    return list_of_data

if __name__ == '__main__':
    lst = loop()
    print(f'found {len(lst)} unique data points')
    np.savetxt('EOF.csv', lst, delimiter=',', header='folderidx,float,profileidx,year,month,day,'
               +'lat,lon,nlat,nlon,ts1,ts2,ts3,ts4,ts5,ts6,ts7,ts8,'
               + 'ssh_slope,ssh_dir,bath_slope,bath_dir,roughness', comments='')
