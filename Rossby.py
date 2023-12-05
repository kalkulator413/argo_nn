from __future__ import print_function
import numpy as np 
from GeneralUtilities.Compute.list import LatList,LonList
from GeneralUtilities.Compute.constants import degree_dist
from GeneralUtilities.Data.Filepath.instance import get_data_folder as get_base_folder
import geopy

class RossbyBase(object):
	def __init__(self):
		dat_file = get_base_folder()+'/Raw/RossbyRad/rossrad.dat'
		data = open(dat_file,'rb')
		lat_list = []
		lon_list = []
		rossby_list = []
		for line in data.readlines():
			lat, lon, rossby_speed, rossby_def = line[:-2].split()
			lat = np.floor(float(lat))
			if np.floor(float(lon))>=180.:
				lon = np.floor(float(lon)-360.)
			else:
				lon = np.floor(float(lon))
			rossby_def = float(rossby_def) # in km
			
			lat_list.append(lat)
			lon_list.append(lon)
			rossby_list.append(rossby_def)
		self.lats = LatList(np.unique(lat_list))
		self.lons = LonList(np.unique(lon_list))
		self.rossby_dict = dict(zip(list(zip(lat_list,lon_list)),rossby_list))

	def return_rossby_def(self,point):
		closest_lat = self.lats.find_nearest(point.latitude)
		closest_lon = self.lons.find_nearest(point.longitude)
		return self.rossby_dict[(closest_lat,closest_lon)]

