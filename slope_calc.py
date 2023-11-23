# inspired from gaussian fit routine found here (https://stackoverflow.com/questions/21566379/fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m)

import scipy.optimize as opt
import numpy as np
from GeneralUtilities.Compute.Depth.depth_utilities import ETopo1Depth
from GeneralUtilities.Compute.constants import calculate_barotropic_rossby_def_rad
from GeneralUtilities.Compute.constants import degree_dist
import geopy
import matplotlib.pyplot as plt


class SlopeCalcDepth(ETopo1Depth):
	def rossby_def_extent(self,point): #uses barotropic deformation radius 
		dist = calculate_barotropic_rossby_def_rad(point)
		ydist = dist/degree_dist
		xdist = dist/(np.cos(np.deg2rad(point.latitude)))
		if xdist > 2:
			xdist = 2
		if ydist >2:
			ydist = 2
		return self.regional_subsample(x0+rdef,x0-rdef,y0+rdef,y0-rdef)

def get_opt_params(subsampled_depth_class):

	def twoD_Gaussian(xy, amplitude, sigma_x, sigma_y, theta, offset):
	    x, y = xy
	    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
	    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
	    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
	    g = offset + amplitude*np.exp( - (a*((x)**2) + 2*b*(x)*(y)
	                            + c*((y)**2)))
	    return g.ravel()

	initial_guess = (3,20,40,0,10) # could and should be optimized, recommend calculating ~ 1000 random locations and choose mean of those 
	x,y = np.meshgrid(subsamp_depth.lon-np.mean(subsamp_depth.lon),subsamp_depth.lat-np.mean(subsamp_depth.lat))
	popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), subsamp_depth.z.data.ravel()-subsamp_depth.z.data.mean(), p0=initial_guess,maxfev=100000) 
	print('The angle is '+str(np.rad2deg(popt[3])%360))
	print(popt)
	data_fitted = twoD_Gaussian((x, y), *popt)
	fig, ax = plt.subplots(1, 1)
	ax.imshow(subsamp_depth.z, cmap=plt.cm.jet, origin='lower',extent=(x.min(), x.max(), y.min(), y.max()))
	endy = y.max() * np.sin(popt[3])
	endx = x.max() * np.cos(popt[3])
	ax.plot([0,endx], [0,endy])
	plt.show()


world_depth = ETopo1Depth.load()
subsamp_depth = world_depth.regional_subsample(144,140,31,28) #radius of this should scale with rossby deformation radius
get_opt_params(subsamp_depth)

subsamp_depth = world_depth.regional_subsample(161,156,38,24) #radius of this should scale with rossby deformation radius
get_opt_params(subsamp_depth)

subsamp_depth = world_depth.regional_subsample(-160,-166,54,52) #radius of this should scale with rossby deformation radius
get_opt_params(subsamp_depth)

