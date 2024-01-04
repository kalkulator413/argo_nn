from __future__ import print_function
import numpy as np 
from GeneralUtilities.Compute.list import LatList,LonList
from GeneralUtilities.Compute.constants import degree_dist
from GeneralUtilities.Data.Filepath.instance import get_data_folder as get_base_folder
from GeneralUtilities.Compute.constants import degree_dist

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
		km_dist = self.rossby_dict[(closest_lat,closest_lon)] # rossby dict return in km
		ydist = km_dist/degree_dist # convert to latitude degree distance
		xdist = km_dist/(degree_dist*np.cos(np.deg2rad(point.latitude))) # convert to longitude degree distance
		return (ydist,xdist)	

	def rossby_def_extent(self,point,grabber): #uses baroclinic deformation radius 
		ydist,xdist = self.return_rossby_def(point)
		if xdist > 2: # if the degree distance is over 2 degrees clip the field
			xdist = 2
		if ydist >2:
			ydist = 2
		return grabber.get_rect(lat=point.latitude, lon=point.longitude, lat_radius=2*ydist, lon_radius=2*xdist)

