import numpy as np
import scipy.io
from bpp import bpp

print('Loading data...')
data = scipy.io.loadmat('irmdata/epinions.mat', squeeze_me=True)
X = np.array(data['X'], dtype='int8')

print('Starting demo...')
rect_locs, row_locs, col_locs, perps = bpp.bpp_mcmc(X)
