import numpy as np
import scipy.sparse.linalg as LA
from argo_nn.Grabber import *

def grid_pca(grid):

    eigenvals, eigenvecs = LA.eigs(grid, k=2)
    # already sorted by magnitude of eigenvalues
    first_pc = eigenvecs[:, 0]
    second_pc = eigenvecs[:, 1]
    print(first_pc.shape)
    # angle?
    
    # each vector has 361 dimensions
    return eigenvecs.shape



test = EtopoGrabber()
grid = test.get_grid(-10, 20, 3)
print(grid_pca(grid))