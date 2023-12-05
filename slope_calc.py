# inspired from gaussian fit routine found here (https://stackoverflow.com/questions/21566379/fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m)

import scipy.optimize as opt
import numpy as np
from GeneralUtilities.Compute.Depth.depth_utilities import ETopo1Depth
from GeneralUtilities.Compute.constants import calculate_barotropic_rossby_def_rad
from GeneralUtilities.Compute.constants import degree_dist
import geopy
import matplotlib.pyplot as plt
from argo_nn.Rossby import RossbyBase

class SlopeCalcDepth(ETopo1Depth):
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
		return self.regional_subsample(point.longitude+xdist,point.longitude-xdist,point.latitude+ydist,point.latitude-ydist)

def get_opt_params(subsample_depth,point):

	x,y = np.meshgrid(subsample_depth.lon-np.mean(subsample_depth.lon),subsample_depth.lat-np.mean(subsample_depth.lat))
	z = subsample_depth.return_z(point)
	H = np.array(list(zip(x.ravel(),y.ravel())))		
	R_inv = np.diag(np.exp(-(x.ravel()**2+y.ravel()**2)))
	P_inv = np.diag([1/5]*2)
	denom = np.linalg.inv(H.T.dot(R_inv.dot(H))+P_inv)
	gain = denom.dot(H.T.dot(R_inv))
	a = gain.dot(subsample_depth.z.ravel()-z)

	fig, ax = plt.subplots(1, 1)
	ax.imshow(subsample_depth.z, cmap=plt.cm.jet, origin='lower',extent=(x.min(), x.max(), y.min(), y.max()))
	ax.plot([0,a[0]/max(a)], [0,a[1]]/max(a))
	plt.show()


world_depth = SlopeCalcDepth.load()
point = geopy.Point(29,142)
subsample_depth = world_depth.rossby_def_extent(point) #radius of this should scale with rossby deformation radius
get_opt_params(subsample_depth,geopy.Point(29,142))

point = geopy.Point(31,159)
subsample_depth = world_depth.rossby_def_extent(point) #radius of this should scale with rossby deformation radius
get_opt_params(subsample_depth,point)

point = geopy.Point(53,-163)
subsample_depth = world_depth.rossby_def_extent(point) #radius of this should scale with rossby deformation radius
get_opt_params(subsample_depth,point)

