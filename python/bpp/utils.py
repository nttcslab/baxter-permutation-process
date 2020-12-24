# utility functions

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import patches
from dataclasses import dataclass

@dataclass(frozen=True)
class BPPOptions:
    """
    Options for Baxter Permutation Process

    Parameters
    ----------
    enc: float
        (to be checked)
    n: int
        Initial length of Baxter permutation
    alpha: float
        (to be checked)
    maxiter: int
        Maximum number of iterations
    missing_ratio: float
        Ratio of missing elements in test data
    rand_seed: int or None
        Fixed random seed, default = None
    """
    enc: float = None
    n: int = 1
    alpha: float = 0.1
    maxiter: int = 1001
    missing_ratio: float = 0.1
    rand_seed: int = None

def extract_maxima(a):
    """
    Searching left-to-right and right-to-left maxima from a given Baxter permutation.

    Parameters
    ----------
    a: array-like
        A Baxter permutation

    Returns
    ----------
    l2r_maxs: numpy.ndarray
        A sequence of left-to-right maxima
    r2l_maxs: numpy.ndarray
        A sequence of right-to-left maxima
    l2r_locs: numpy.ndarray
        A sequence of indices of the left-to-right maxima
    r2l_locs: numpy.ndarray
        A sequence of indices of the right-to-left maxima
    """

    # initialization
    n = a.shape[0]
    l2r_maxs = [ a[0] ]
    r2l_maxs = [ a[-1] ]
    l2r_locs = [ 0 ]
    r2l_locs = [ a.shape[0]-1 ]

    # search
    for i_left in range(1, n):
        i_right = (n-1)-i_left
        ## left
        if a[i_left] > l2r_maxs[-1]:
            l2r_maxs.append(a[i_left])
            l2r_locs.append(i_left)
        # right
        if a[i_right] > r2l_maxs[-1]:
            r2l_maxs.append(a[i_right])
            r2l_locs.append(i_right)

    # return
    l2r_maxs = np.array(l2r_maxs)
    r2l_maxs = np.array(r2l_maxs)
    l2r_locs = np.array(l2r_locs)
    r2l_locs = np.array(r2l_locs)
    return l2r_maxs, r2l_maxs, l2r_locs, r2l_locs

def plot_state(figs, X, rect_locs, row_locs, col_locs, perps):
    """
    Visualizing the current result

    Parameters
    ----------
    figs: tuple
        A tuple of a pyplot Figure object and pyploy Axes objects
    X: numpy.ndarray
        Input data
    rect_locs: numpy.ndarray
        A matrix representing a rectangular partitioning
    row_locs: numpy.ndarray
        A vector representing row locations of rectangles in input data
    col_locs: numpy.ndarray
        A vector representing column locations of rectangles in input data
    perps: numpy.ndarray
        A sequence of test perplexity values
    """

    # initializing
    fig, (ax1, ax2) = figs
    height, width = X.shape

    # Sorting
    sorted_row_indices = np.argsort(row_locs)
    sorted_col_indices = np.argsort(col_locs)
    sorted_row_locs = row_locs[ sorted_row_indices ]
    sorted_col_locs = col_locs[ sorted_col_indices ]

    # showing sorted input data
    ax1.clear()
    tmp_X = X[sorted_row_indices]
    tmp_X = tmp_X[:, sorted_col_indices]
    ax1.matshow(tmp_X, cmap='binary')

    # drawing perplexity
    ax2.clear()
#    ax2.set_ylim([0, 5])
    ax2.set_yscale('log')
    ax2.set_xlabel('MCMC iterations')
    ax2.set_ylabel('Perplexity')
    ax2.plot(perps)

    # drawing rectangles
    num_blocks = rect_locs.shape[0]
    for ii in range(num_blocks):
        cur_rect_loc = rect_locs[ii]
        ##
        cnd_r = (cur_rect_loc[0] < sorted_row_locs) & (sorted_row_locs <= cur_rect_loc[1])
        extracted_rows = sorted_row_locs[cnd_r]
        ##
        cnd_c = (cur_rect_loc[2] < sorted_col_locs) & (sorted_col_locs <= cur_rect_loc[3])
        extracted_cols = sorted_col_locs[cnd_c]
        ##
        if (len(extracted_rows)>0) and (len(extracted_cols)>0):
            ##
            min_r, max_r = extracted_rows.min()*height-0.5, extracted_rows.max()*height+0.5
            min_c, max_c = extracted_cols.min()*width -0.5, extracted_cols.max()*width +0.5
            rect_width, rect_height = max_c-min_c, max_r-min_r
            ## drawing the current rectangle
            draw_rect = patches.Rectangle(xy=(min_c, min_r),
                                          width=rect_width, height=rect_height,
                                          fill=False, color='blue', linewidth=3)
            ax1.add_patch(draw_rect)

    # drawing
    plt.draw(); plt.pause(0.01)
