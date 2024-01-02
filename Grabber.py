from GeneralUtilities.Data.Filepath.instance import get_data_folder as get_base_folder
from netCDF4 import Dataset
import numpy as np
from math import isnan

class GridHolder():
	def __init__(self, latlist, lonlist, grid):
		self.latlist = latlist
		self.lonlist = lonlist
		self.grid = grid

class Grabber():

	lat =  None
	lon = None

	def get_rect(self, lat, lon, lat_radius, lon_radius):
		### exclude edge cases near +- 90, +- 180
		assert lat - lat_radius > -90 and lat + lat_radius < 90, 'latitude out of bounds'
		assert lon - lon_radius > -180 and lon + lon_radius < 180, 'longtitude out of bounds'

		lon_dist = round(lon_radius / self.granularity)
		lat_dist = round(lat_radius / self.granularity)
		lat_idx = round((lat + 90 - self.offset) / self.granularity)
		lon_idx = round((lon + 180 - self.offset) / self.granularity)
		
		result = self.z[lat_idx - lat_dist:lat_idx + lat_dist + 1, lon_idx - lon_dist:lon_idx + lon_dist + 1]
		### exclude coastal areas (Aviso)
		assert not result.mask.any(), 'tried to grab ssh over land'

		latlist = self.lat[lat_idx - lat_dist:lat_idx + lat_dist + 1]
		lonlist = self.lon[lon_idx - lon_dist:lon_idx + lon_dist + 1]
		return GridHolder(latlist, lonlist, (result).data)

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
