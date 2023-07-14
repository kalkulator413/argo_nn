from GeneralUtilities.Data.Filepath.instance import get_data_folder as get_base_folder
from netCDF4 import Dataset
import numpy as np
from math import isnan

class Grabber():
	def get_grid(self, lat:float, lon:float, radius:float):

		### exclude edge cases near +- 90, +- 180
		assert lat - radius > -90 and lat + radius < 90, 'latitude out of bounds'
		assert lon - radius > -180 and lat + radius < 180, 'longtitude out of bounds'

		dist = round(radius / self.granularity)
		lat_idx = round((lat + 90 - self.offset) / self.granularity)
		lon_idx = round((lon + 180 - self.offset) / self.granularity)
		
		result = np.zeros( (dist * 2 + 1) * (dist * 2 + 1))
		ctr = 0
		for lt in range(lat_idx - dist, lat_idx + dist + 1):
			for ln in range(lon_idx - dist, lon_idx + dist + 1):
				result[ctr] = self.z[lt][ln]
				### exclude coastal areas (Aviso)
				assert not isnan(result[ctr]), 'tried to grab ssh over land'

				ctr += 1
		
		return result.reshape((dist * 2 + 1, dist * 2 + 1))

class EtopoGrabber(Grabber):
	def __init__(self):
		nc_fid_z_data = Dataset(get_base_folder()+'/Raw/ETopo1/ETOPO1_Bed_c_gdal.grd')
		nc_fid_coord = Dataset(get_base_folder()+'/Raw/ETopo1/ETOPO1_Bed_g_gmt4.grd')
		self.lon = nc_fid_coord['x'][:-1].data
		self.lat = nc_fid_coord['y'][:-1].data
		self.z = nc_fid_z_data['z'][:].reshape(len(self.lat),len(self.lon))
		self.granularity = 0.01666666666666
		self.offset = 0

class AvisoGrabber(Grabber):
	def __init__(self, year:int, month:int, day:int):

		# date format checks
		assert day > 0 and day <= 31, 'day must be between 1 and 31'
		assert month > 0 and month <= 12, 'month must be between 1 and 12'
		assert year >= 2017 and year <= 2023, 'year must be between 2017 and 2023'
		if year == 2017:
			assert month > 2, 'jan, feb 2017 not in data set'
		elif year == 2021:
			assert month != 12, 'dec 2021 not in data set'
		elif year == 2022:
			assert month != 1, 'jan 2022 not in data set'
		elif year == 2023:
			assert month < 7, 'this month in 2023 has not occured yet'

		nc_ssh = Dataset(get_base_folder() + f'/Raw/Aviso/{year}/{month}.nc')
		self.lat = nc_ssh['latitude'][:].data
		self.lon = nc_ssh['longitude'][:].data
		self.z = nc_ssh['sla'][day-1]
		self.granularity = 0.25
		self.offset = 0.125
