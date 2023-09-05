from GeneralUtilities.Data.Filepath.instance import get_data_folder as get_base_folder
from netCDF4 import Dataset
import numpy as np
import pickle

sal_file = get_base_folder() + '/Raw/Argo/temp-sal/RG_ArgoClim_Salinity_2019.nc'
temp_file = get_base_folder() + '/Raw/Argo/temp-sal/RG_ArgoClim_Temperature_2019.nc'
sal_ds = Dataset(sal_file)
temp_ds = Dataset(temp_file)
pressure_range = range(25, 56)
mean_sal = [sal_ds['ARGO_SALINITY_MEAN'][depth] for depth in pressure_range]
mean_temp = [temp_ds['ARGO_TEMPERATURE_MEAN'][depth] for depth in pressure_range]

class TS_Grabber():
    def __init__(self, year:int, month:int):
        assert month > 0 and month <= 12, 'month must be between 1 and 12'

        if year <= 2018:
            # calculate how many months back it is from dec 2018
            month_diff = (year*12 + month) - (2018*12 + 12) - 1

            self.temp = [temp_ds['ARGO_TEMPERATURE_ANOMALY'][month_diff][depth] for depth in pressure_range]
            self.sal = [sal_ds['ARGO_SALINITY_ANOMALY'][month_diff][depth] for depth in pressure_range]
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

            self.temp = [ds['ARGO_TEMPERATURE_ANOMALY'][0][depth] for depth in pressure_range]
            self.sal = [ds['ARGO_SALINITY_ANOMALY'][0][depth] for depth in pressure_range]

        self.lat = temp_ds['LATITUDE'][:]
        self.lon = temp_ds['LONGITUDE'][:]

    def get_profiles(self, lat:float, lon:float):
        assert lat >= -64.5 and lat < 80, 'given latitude not in data'

        if lon < 20:
            lon += 360
        
        lat_idx = round(lat + 64.5)
        lon_idx = round(lon - 20.5)

        temps = np.zeros(len(pressure_range), dtype='float32')
        sals = np.zeros(len(pressure_range), dtype='float32')
        for i in range(len(self.temp)):

            ta = self.temp[i].data[lat_idx][lon_idx]
            sa = self.sal[i].data[lat_idx][lon_idx]
            tm = mean_temp[i].data[lat_idx][lon_idx]
            sm = mean_sal[i].data[lat_idx][lon_idx]

            assert ta != -999 and sa != -999 and tm != -999 and sm != -999, 'used masked data'
            temps[i] = ta + tm
            sals[i] = sa + sm

        return temps, sals

def get_eofs():
    time_len = 234
    N = 397 * time_len
    L = len(pressure_range) * 2
    Y = np.zeros((N, L), dtype='float32')

    sum_time = 0
    for year in range(2004, 2024):
        for month in range(1, 13):
            if (year == 2023 and month > 6):
                break
            
            data = TS_Grabber(year, month)

            sum_coord = 0
            for lat in data.lat[::9]:
                for lon in data.lon[::9]:
                    try:
                        t, s = data.get_profiles(lat, lon)
                    except Exception as e:
                        continue

                    Y[sum_coord * time_len + sum_time] = np.append(t, s)
                    sum_coord += 1

            sum_time += 1

    Y = Y.T
    means = np.zeros(2 * len(pressure_range))
    for i, y in enumerate(Y):
        means[i] = np.mean(y)
        Y[i] = y - means[i]

    with open('eofs/means.pkl', 'wb') as f:
        pickle.dump(means, f)

    print('Y has dimensions', Y.shape)
    U, amp, var = get_decomp(Y)

    with open('eofs/u.pkl', 'wb') as f:
        pickle.dump(U, f)

    with open('eofs/amp.pkl', 'wb') as f:
        pickle.dump(amp, f)

    with open('eofs/var.pkl', 'wb') as f:
        pickle.dump(var, f)

def get_decomp(Y):
    U, s, Vh = np.linalg.svd(Y)

    print(U.shape, s.shape, Vh.shape)
    amp = U.T @ Y
    variances = (1/Y.shape[1]) * s * s

    print('matrix of EOFs has dimensions', U.shape)
    print('matrix of amplitudes has dimensions', amp.shape)
    print('vector of variances has dimensions', variances.shape)

    return U, amp, variances

if __name__ == '__main__':
    get_eofs()