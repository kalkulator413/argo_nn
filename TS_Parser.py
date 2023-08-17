from GeneralUtilities.Data.Filepath.instance import get_data_folder as get_base_folder
from netCDF4 import Dataset
import numpy as np
from math import isnan
import pickle

# find means
sal_file = get_base_folder() + '/Raw/Argo/temp-sal/RG_ArgoClim_Salinity_2019.nc'
temp_file = get_base_folder() + '/Raw/Argo/temp-sal/RG_ArgoClim_Temperature_2019.nc'
sal_ds = Dataset(sal_file)
mean_sal = [sal_ds['ARGO_SALINITY_MEAN'][depth] for depth in range(58)]
temp_ds = Dataset(temp_file)
mean_temp = [temp_ds['ARGO_TEMPERATURE_MEAN'][depth] for depth in range(58)]

class TS_Grabber():
    def __init__(self, year:int, month:int):
        assert month > 0 and month <= 12, 'month must be between 1 and 12'

        if year <= 2018:
            # calculate how many months back it is from dec 2018
            month_diff = (year*12 + month) - (2018*12 + 12) - 1

            self.temp = [temp_ds['ARGO_TEMPERATURE_ANOMALY'][month_diff][depth] for depth in range(58)]
            self.sal = [sal_ds['ARGO_SALINITY_ANOMALY'][month_diff][depth] for depth in range(58)]
        else:
            assert year <= 2023, 'provided year is in the future'
            if year == 2023:
                assert month <= 6, 'this month in 2023 has not happened yet'
            
            path = get_base_folder() + '/Raw/Argo/temp-sal/'
            if month >= 10:
                path += f'RG_ArgoClim_{year}{month}_2019.nc'
            else:
                path += f'RG_ArgoClim_{year}0{month}_2019.nc'
            ds = Dataset(path)

            self.temp = [ds['ARGO_TEMPERATURE_ANOMALY'][0][depth] for depth in range(58)]
            self.sal = [ds['ARGO_SALINITY_ANOMALY'][0][depth] for depth in range(58)]

        self.lat = temp_ds['LATITUDE'][:]
        self.lon = temp_ds['LONGITUDE'][:]

    def get_profiles(self, lat:float, lon:float):
        assert lat >= -64.5 and lat < 80, 'given latitude not in data'

        if lon < 20:
            lon += 360
        
        lat_idx = round(lat + 64.5)
        lon_idx = round(lon - 20.5)

        temps = np.zeros(58)
        sals = np.zeros(58)
        for i in range(58):
            temps[i] = self.temp[i].data[lat_idx][lon_idx] + mean_temp[i].data[lat_idx][lon_idx]
            sals[i] = self.sal[i].data[lat_idx][lon_idx]  + mean_sal[i].data[lat_idx][lon_idx]

            assert not self.temp[i].mask[lat_idx][lon_idx], 'no data'
            assert not self.sal[i].mask[lat_idx][lon_idx], 'no data'
            assert not mean_temp[i].mask[lat_idx][lon_idx], 'no data'
            assert not mean_sal[i].mask[lat_idx][lon_idx], 'no data'

        return temps, sals

def get_eofs():
    N = 48160
    L = 58
    Y_t = np.zeros((N, L))
    Y_s = np.zeros((N, L))

    sum = 0
    for year in range(2004, 2024):
        for month in range(1, 12):
            if (year == 2023 and month > 6):
                break
            
            data = TS_Grabber(year, month)

            for lat in data.lat[::12]:
                for lon in data.lon[::12]:
                    try:
                        t, s = data.get_profiles(lat, lon)
                    except Exception as e:
                        continue

                    Y_t[sum] = t
                    Y_s[sum] = s
                    sum += 1

    print('sum:', sum)
    print('temperature matrix has dimensions', Y_t.shape)
    U_t, amp_t, var_t = get_decomp(Y_t)
    
    print('\nsalinity matrix has dimensions', Y_s.shape)
    U_s, amp_s, var_s = get_decomp(Y_s)

    with open('testing-notebooks/ut_full.pkl', 'wb') as f:
        pickle.dump(U_t, f)

    with open('testing-notebooks/us_full.pkl', 'wb') as f:
        pickle.dump(U_s, f)

def get_decomp(Y):
    U, s, Vh = np.linalg.svd(Y)
    amp = np.diag(s) @ Vh
    variances = np.diag((np.diag(s) @ np.diag(s).T)) / 58

    print('matrix of EOFs has dimensions', U.shape)
    print('matrix of amplitudes has dimensions', amp.shape)
    print('vector of variances has dimensions', variances.shape)
    print('now printing variances')
    for v in variances:
        print(v)

    return U, amp, variances

# get_eofs()