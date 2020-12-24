import scipy.io
from bpp import bpp

print('Loading data...')
data = scipy.io.loadmat('irmdata/epinions.mat', squeeze_me=True)
X = data['X']

print('Starting demo...')
rect_locs, row_locs, col_locs, perps = bpp.bpp_mcmc(X)
