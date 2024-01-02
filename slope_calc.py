# inspired from gaussian fit routine found here (https://stackoverflow.com/questions/21566379/fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m)

import numpy as np
from GeneralUtilities.Compute.constants import degree_dist
import geopy
import matplotlib.pyplot as plt
from argo_nn.Rossby import RossbyBase
from argo_nn.Grabber import *

class SlopeCalcDepth(EtopoGrabber):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)		
		self.rossby = RossbyBase()

	def rossby_def_extent(self,point): #uses baroclinic deformation radius 
		dist = self.rossby.return_rossby_def(point)
		ydist = dist/degree_dist
		xdist = dist/(np.cos(np.deg2rad(point.latitude)))
		if xdist > 2:
			xdist = 2
		if ydist >2:
			ydist = 2
		return self.get_rect(point.latitude, point.longitude, lat_radius=ydist, lon_radius=xdist)
	
class SlopeCalcSSH(AvisoGrabber):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)		
		self.rossby = RossbyBase()

	def rossby_def_extent(self,point): #uses baroclinic deformation radius 
		dist = self.rossby.return_rossby_def(point)
		ydist = dist/degree_dist
		xdist = dist/(np.cos(np.deg2rad(point.latitude)))
		if xdist > 2:
			xdist = 2
		if ydist >2:
			ydist = 2
		return self.get_rect(point.latitude, point.longitude, lat_radius=ydist, lon_radius=xdist)

# TODO: scale coviariance matrices by Rossby deformation radius

def get_opt_params(subsample_depth,point):

	x,y = np.meshgrid(subsample_depth.lonlist-np.mean(subsample_depth.lonlist),subsample_depth.latlist-np.mean(subsample_depth.latlist))
	z = subsample_depth.grid[len(subsample_depth.latlist) // 2][len(subsample_depth.lonlist) // 2]
	H = np.array(list(zip(x.ravel(),y.ravel())))		
	R_inv = np.diag(np.exp(-(x.ravel()**2+y.ravel()**2)))
	P_inv = np.diag([1/5]*2)
	denom = np.linalg.inv(H.T.dot(R_inv.dot(H))+P_inv)
	gain = denom.dot(H.T.dot(R_inv))
	a = gain.dot(subsample_depth.grid.ravel()-z)

	fig, ax = plt.subplots(1, 1)
	ax.imshow(subsample_depth.grid, cmap=plt.cm.jet, origin='lower',extent=(x.min(), x.max(), y.min(), y.max()))
	ax.plot([0,a[0]/max(a)], [0,a[1]]/max(a))
	plt.show()

if __name__ == '__main__':
	world_depth = SlopeCalcSSH(2017, 11, 8)
	# world_depth = SlopeCalcDepth()
	point = geopy.Point(29,142)
	subsample_depth = world_depth.rossby_def_extent(point) #radius of this should scale with rossby deformation radius
	get_opt_params(subsample_depth,point)

	point = geopy.Point(31,159)
	subsample_depth = world_depth.rossby_def_extent(point) #radius of this should scale with rossby deformation radius
	get_opt_params(subsample_depth,point)

	point = geopy.Point(53,-163)
	subsample_depth = world_depth.rossby_def_extent(point) #radius of this should scale with rossby deformation radius
	get_opt_params(subsample_depth,point)

