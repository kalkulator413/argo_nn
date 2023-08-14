from GeneralUtilities.Data.Filepath.instance import get_data_folder as get_base_folder
from netCDF4 import Dataset
import numpy as np
from math import isnan

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
    # N = ((2022 - 2004) * 12 + 6) * 30650 # 30650 is the total number of valid pairs (lat, lon) with unmasked data
    N = 30650
    L = 58
    temp_Y = np.zeros((N, L))
    sal_Y = np.zeros((N, L))

    sum = 0
    for year in [2004]:#range(2004, 2024):
        for month in [1]:#range(1, 12):
            if (year == 2023 and month > 6):
                break
            
            data = TS_Grabber(year, month)

            for lat in data.lat:
                for lon in data.lon:
                    try:
                        t, s = data.get_profiles(lat, lon)
                    except Exception as e:
                        continue

                    temp_Y[sum] = t
                    sal_Y[sum] = s
                    sum += 1

    print(temp_Y.shape)
    U_t, S_t, V_t = np.linalg.svd(temp_Y)
    amp_t = S_t @ V_t.T
    print(U_t.shape, amp_t.shape)

    # U_s, S_s, V_s = np.linalg.svd(sal_Y)
    # amp_s = U_s.T @ sal_Y

get_eofs()