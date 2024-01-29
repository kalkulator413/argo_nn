# inspired from gaussian fit routine found here (https://stackoverflow.com/questions/21566379/fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m)

import numpy as np
import geopy
import matplotlib.pyplot as plt
from argo_nn.Rossby import RossbyBase
from argo_nn.Grabber import *

# TODO: scale coviariance matrices by Rossby deformation radius

def get_opt_params(subsample_depth, rossby, name, plot=True):

	rossby_xdist, rossby_ydist = rossby.return_rossby_def(geopy.Point(np.mean(subsample_depth.latlist),np.mean(subsample_depth.lonlist)))	
	x,y = np.meshgrid(subsample_depth.lonlist-np.mean(subsample_depth.lonlist),subsample_depth.latlist-np.mean(subsample_depth.latlist))
	z = subsample_depth.grid[len(subsample_depth.latlist) // 2][len(subsample_depth.lonlist) // 2]

	H = np.array(list(zip(x.ravel(),y.ravel())))		
	R_inv = np.diag(np.exp(-((x.ravel()/rossby_xdist)**2+(y.ravel()/rossby_ydist)**2))) # R sets the observational noise
	P_inv = np.diag(np.array([1, 1]))**2 # ratio of P to R is noise to signal
	denom = np.linalg.inv(H.T @ R_inv @ H + P_inv)
	gain = R_inv @ H @ denom
	a = (subsample_depth.grid.ravel()-z) @ gain

	if plot:
		fig, ax = plt.subplots(1, 1)
		ax.imshow(subsample_depth.grid, cmap=plt.cm.jet, origin='lower',extent=(x.min(), x.max(), y.min(), y.max()))
		ax.plot([0,a[0]/max(np.abs(a))], [0,a[1]/max(np.abs(a))])
		plt.show()

		# for when i'm using the remote server, because it doesnt show the plots normally
		# plt.savefig(name)

	return a[0]/max(np.abs(a)), a[1]/max(np.abs(a))

if __name__ == '__main__':
	world_depth = AvisoGrabber(2017, 11, 8)
	# world_depth = EtopoGrabber()
	rossby = RossbyBase()

	point = geopy.Point(29,142)
	subsample_depth = rossby.rossby_def_extent(point, world_depth) #radius of this should scale with rossby deformation radius
	get_opt_params(subsample_depth, rossby, '1.png')

	point = geopy.Point(31,159)
	subsample_depth = rossby.rossby_def_extent(point, world_depth) #radius of this should scale with rossby deformation radius
	get_opt_params(subsample_depth, rossby, '2.png')

	point = geopy.Point(52,-163)
	subsample_depth = rossby.rossby_def_extent(point, world_depth) #radius of this should scale with rossby deformation radius
	get_opt_params(subsample_depth, rossby, '3.png')

	point = geopy.Point(25.431519, -89.523394)
	subsample_depth = rossby.rossby_def_extent(point, world_depth) #radius of this should scale with rossby deformation radius
	get_opt_params(subsample_depth, rossby, '4.png')

