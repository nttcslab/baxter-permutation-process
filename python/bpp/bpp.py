# Baxter permutation process

import numpy as np
from numpy import random
import scipy.special
import matplotlib.pyplot as plt
import math
import sys
from . import utils

def delete_max(bp, uniform_random):
    """
    Deleting the maximum element from a given Baxter permutation

    Parameters
    ----------
    bp: numpy.ndarray
        A Baxter permutation
    uniform_random: numpy.ndarray
        A sequence of random variables, with the same size as "bp".

    Returns
    ----------
    new_bp: numpy.ndarray
        An updated Baxter permutation
    new_uniform_random:
        An updated random variable sequence.

    """

    # max element
    i = bp.argmax()
    a = bp[i]
    # deleting the max element
    new_bp = np.delete(bp, i)
    new_uniform_random = uniform_random[:-1].copy()

    return new_bp, new_uniform_random

def evolv_bp(bp, uniform_rand, new_uniform_rand=None):
    """
    Adding a new index to a given Baxter permutation

    Parameters
    ----------
    bp: numpy.ndarray
        A Baxter permutation
    uniform_rand: numpy.ndarray
        A sequence of random variables, with the same size as "bp".
    new_uniform_rand: float
        A random variable

    Returns
    ----------
    new_bp: numpy.ndarray
        An updated "bp"
    uniform_rand:
        An updated "uniform_rand"

    Notes
    ----------
    1. If 'new_uniform_rand' is not given or is None, 'new_uniform_rand' will be generated inside,
       and append it to the end of 'uniform_rand'.
    2. If 'new_uniform_rand' is given, we suppose that you have already had a full length of
       'uniform_rand', and thus 'new_uniform_rand' will be used as-is.
    """

    # initialization
    len_bp = bp.shape[0]
    if new_uniform_rand is None:
        new_uniform_rand = random.uniform(0, 1)

    # extracting maxima
    l2r_maxs, r2l_maxs, l2r_locs, r2l_locs = utils.extract_maxima(bp)

    # searching for a position to insert a new index
    # Note: Only a single condition will be accepted among the following 4 ones.
    ## first from left
    if new_uniform_rand < uniform_rand[ l2r_maxs[0] ]:
        new_bp = np.insert(bp.copy(), l2r_locs[0], len_bp)
    ## later from left
    elif l2r_maxs.shape[0] > 1:
        for ii in range(1, l2r_maxs.shape[0]):
            if new_uniform_rand < uniform_rand[ l2r_maxs[ii] ]:
                new_bp = np.insert(bp.copy(), l2r_locs[ii], len_bp); break
    ## first from right
    if ('new_bp' not in locals()) and (uniform_rand[ r2l_maxs[0] ] < new_uniform_rand):
        new_bp = np.insert(bp.copy(), r2l_locs[0]+1, len_bp)
    elif r2l_maxs.shape[0] > 1:
        for ii in range(1, r2l_maxs.shape[0]):
            if uniform_rand[ r2l_maxs[ii] ] < new_uniform_rand:
                new_bp = np.insert(bp.copy(), r2l_locs[ii]+1, len_bp); break
    if 'new_bp' not in locals():
        print('Error: Cannot insert a new element in the existing Baxter permutation.'); sys.exit(1)

    # return
    return new_bp, np.append(uniform_rand, new_uniform_rand)

def bp2fp(bp):
    """
    Converting a given Baxter permutation to a floorplan partitioning

    Parameters
    ----------
    bp: numpy.ndarray
        A Baxter permutation

    Returns
    ----------
    label_mat: numpy.ndarray
        A matrix representing a floorplan partitioning corresponding to the Baxter permutation
    """

    # initializing
    num_blocks = bp.shape[0]
    label_mat = np.ones((num_blocks, num_blocks), dtype=int) * bp[0]
    ## blocks
    blocks = np.zeros((num_blocks, 4), dtype=int)  # block representation = [top, bottom, left, right]
    blocks[bp[0]] = [0, num_blocks, 0, num_blocks]

    # for every index in the permutation
    for ii in range(1, num_blocks):
        cur_block = blocks[bp[ii]]
        prev_block = blocks[bp[ii-1]]
        ## if the current index is smaller than the previous one
        ##   --> creating a new block at the right of the previous one
        if bp[ii-1] < bp[ii]:
            ## allocating the current block
            blocks[bp[ii]] = [ prev_block[0], prev_block[1], ii, prev_block[3] ]
            label_mat[ prev_block[0]:prev_block[1], ii:prev_block[3] ] = bp[ii]
            ## updating the previous (= left) block
            blocks[bp[ii-1]] = [ prev_block[0], prev_block[1], prev_block[2], ii ]
            ## checking the under neighbor block
            bottom_cur_block = cur_block[1]  # bottom of the current block
            while bottom_cur_block < num_blocks:
                cur_block = blocks[bp[ii]]
                ## check only if the index of the under block is smaller than the one of the current block
                under_block_index = label_mat[ bottom_cur_block, cur_block[2] ]  # index of the under neightbor
                under_block = blocks[under_block_index]  # the under neightbor block
                if under_block_index > bp[ii]: break
                ## extending the current block below, to the bottom of the under neighbor block
                blocks[bp[ii], 1] = under_block[1]
                label_mat[ cur_block[0]:under_block[1], cur_block[2]:cur_block[3] ] = bp[ii]
                bottom_cur_block = under_block[1]
                ## shrinking the under neighbor block to left, to the left edge of the current block
                blocks[under_block_index, 3] = cur_block[2]
        ## else --> creating a new block above the previous one
        else:
            ## allocating the current block
            blocks[bp[ii]] = [ prev_block[0], num_blocks-ii, prev_block[2], prev_block[3] ]
            label_mat[ prev_block[0]:num_blocks-ii, prev_block[2]:prev_block[3] ] = bp[ii]
            ## updating the previous (= the under neighbor) block
            blocks[bp[ii-1]] = [ num_blocks-ii, prev_block[1], prev_block[2], prev_block[3] ]
            ## checking the left block
            left_cur_block = cur_block[2]
            while left_cur_block > 0:
                cur_block = blocks[bp[ii]]
                ## check only if the index of the left block is larger than the of the current block
                left_block_index = label_mat[cur_block[0], left_cur_block-1]
                left_block = blocks[left_block_index]
                if left_block_index < bp[ii]: break
                ## extending the current block to the left edge of the left block
                blocks[bp[ii], 2] = left_block[2]
                label_mat[ cur_block[0]:cur_block[1], left_block[2]:cur_block[3] ] = bp[ii]
                left_cur_block = left_block[2]
                ## shrinking the left block below, to the bottom of the current block
                blocks[left_block_index, 0] = cur_block[1]

    return label_mat

def labelmat2blocks(label_mat):
    """
    Converting a given label matrix to a set of blocks

    Parameters
    ----------
    label_matr: numpy.ndarray
        A matrix representing a floorplan partitioning

    Returns
    ----------
    blocks: numpy.ndarray
        A set of blocks in a given label matrix
    """

    num_blocks = label_mat.max()+1
    blocks = np.zeros((num_blocks, 4), dtype=int)
    for ii in range(num_blocks):
        check = np.where(label_mat==ii)
        blocks[ii] = [ check[0].min(), check[0].max()+1, check[1].min(), check[1].max()+1 ]
    return blocks

def evolv_rp(bp, label_mat, beta_rand, rect_locs):
    """
    Adding a new rectangle according to a given "label matrix"

    Parameters
    ----------
    bp: numpy.ndarray
        A Baxter permutation
    label_mat: numpy.ndarray
        A matrix representing a floorplan partitioning corresponding to the Baxter permutation
    beta_rand: float
        A random variable drawn from a Beta distribution
    rect_locs: numpy.ndarray
        A matrix representing a rectangular partitioning

    Returns
    ----------
    new_rect_locs: numpy.ndarray
        updated "rect_locs"
    """

    # initialization
    ## a set of blocks
    blocks = labelmat2blocks(label_mat)
    ## allocating a new rectangle for the largest index in "label_matrix"
    new_rect_locs = np.append(rect_locs.copy(), np.zeros((1, 4)), axis=0)

    prev_block = blocks[-2]  # a block having the second largest block index
    cur_block = blocks[-1]   # a block having the largest block index
    # if "cur_block" is below "prev_block"
    if prev_block[1] < cur_block[1]:
        # the above blocks and rectangles
        above_blk_idxs = np.unique( label_mat[ cur_block[0]-1, cur_block[2]:cur_block[3] ] )
        above_rects = rect_locs[above_blk_idxs]
        # new height of the current rectangle
        height_above_rect = np.min(above_rects[:,1] - above_rects[:,0])
        cut_length = (1 - beta_rand) * height_above_rect
        # creating a new rectangle below
        new_rect_locs[above_blk_idxs, 1] = above_rects[:,1] - cut_length
        new_rect_locs[-1, 0] = above_rects.max(axis=0)[1] - cut_length
        new_rect_locs[-1, 1] = above_rects.max(axis=0)[1]
        new_rect_locs[-1, 2] = above_rects.min(axis=0)[2]
        new_rect_locs[-1, 3] = above_rects.max(axis=0)[3]
    else:
        # the left block and rectangle
        left_blk_idxs = np.unique( label_mat[ cur_block[0]:cur_block[1], cur_block[2]-1 ] )
        left_rects = rect_locs[left_blk_idxs]
        # new width of the current rectangle
        width_left_rect = np.min(left_rects[:,3] - left_rects[:,2])
        cut_length = (1 - beta_rand) * width_left_rect
        # creating a new rectangle right
        new_rect_locs[left_blk_idxs, 3] = left_rects[:,3] - cut_length
        new_rect_locs[-1, 0] = left_rects.min(axis=0)[0]
        new_rect_locs[-1, 1] = left_rects.max(axis=0)[1]
        new_rect_locs[-1, 2] = left_rects.max(axis=0)[3] - cut_length
        new_rect_locs[-1, 3] = left_rects.max(axis=0)[3]

    return new_rect_locs

def delete_rp(bp, label_mat, rect_locs):
    """
    Deleting a rectangle having the largest index by enlarging a neighboring rectangle

    Parameters
    ----------
    bp: numpy.ndarray
        A Baxter permuatation
    label_mat: numpy.ndarray
        A matrix representing a floorplan partitioning corresponding to the Baxter permutation
    rect_locs: numpy.ndarray
        A matrix representing a rectangular partitioning

    Returns
    ----------
    rect_locs:
        An updated "rect_locs"
    """

    # initialization
    ## a set of blocks
    blocks = labelmat2blocks(label_mat)

    prev_block = blocks[-2]  # a block having the second largest block index
    cur_block = blocks[-1]   # a block having the largest block index
    cur_rect = rect_locs[-1]  # a rectangle having the largest block index
    # if "cur_block" is below "prev_block"
    if prev_block[1] < cur_block[1]:
        # deleting a rectangle by enlarging the above rectangle
        above_blk_idxs = label_mat[ cur_block[0]-1, cur_block[2]:cur_block[3] ]
        height_cur_rect = cur_rect[1] - cur_rect[0]
        rect_locs[above_blk_idxs, 1] += height_cur_rect
    else:
        # deleting a rectangle by enlarging the left rectangle
        left_blk_idxs = label_mat[ cur_block[0]:cur_block[1], cur_block[2]-1 ]
        width_cur_rect = cur_rect[3] - cur_rect[2]
        rect_locs[left_blk_idxs, 3] += width_cur_rect

    return rect_locs

def count_elems_in_rects(X, rect_locs, row_locs, col_locs, row_pos=None, col_pos=None):
    """
    Counting the number of elements having a specific value for each rectangle

    Parameters
    ----------
    X: numpy.ndarray
        Input data
    rect_locs: numpy.ndarray
        Rectangular partitioning
    row_locs: numpy.ndarray
        A vector representing row locations of rectangles in input data
    col_locs: numpy.ndarray
        A vector representing column locations of rectangles in input data
    row_pos: array-like or None
        A row position in the input data to be checked, None means all.
    col_pos: array-like or None
        A column position in the input data to be checked, None means all.

    Returns
    ----------
    each_rect_count: numpy.ndarray
        Number of elements having a value [row] for each rectangle [col]
    """

    # initialization
    ## Suppose that
    ## (1) X has integer elements
    ## (2) missing elements are represented as -1
    max_elems = X.max() + 1
    num_blocks = rect_locs.shape[0]
    height = row_locs.shape[0]
    width  = col_locs.shape[0]
    each_rect_count = np.zeros((max_elems, num_blocks))
    row_pos_bool = np.ones(height, dtype='bool')
    col_pos_bool = np.ones(width,  dtype='bool')
    if row_pos is not None:
        row_pos_bool = np.logical_not(row_pos_bool)
        row_pos_bool[row_pos] = True
    if col_pos is not None:
        col_pos_bool = np.logical_not(col_pos_bool)
        col_pos_bool[col_pos] = True

    # for each block
    for ii in range(num_blocks):
        ## location of the current rectangle
        cur_rect_loc = rect_locs[ii]
        ## searching for a block of input data, which belongs to the current rectangle
        _row_pos = (row_locs > cur_rect_loc[0]) & (row_locs <= cur_rect_loc[1]) & (row_pos_bool)
        _col_pos = (col_locs > cur_rect_loc[2]) & (col_locs <= cur_rect_loc[3]) & (col_pos_bool)
        extracted_data = X[_row_pos][:,_col_pos]
        ## checking the number of elements
        for jj in range(max_elems):
            each_rect_count[jj, ii] = np.count_nonzero(extracted_data==jj)

    return each_rect_count

def calc_Dirichlet_likelihood(each_rect_count, alpha):
    """
    Computing a Dirichlet log likelihood for a given

    Parameters
    ----------
    each_rect_count: numpy.ndarray
        Number of elements having a value [row] for each rectangle [col]
    alpha: float
        A Dirichlet parameter

    Returns
    ----------
    likelihood : float
        Likelihood
    """

    # initialization
    max_elems = each_rect_count.shape[0]
    num_blocks = each_rect_count.shape[1]
    # computing a Likelihood
    ## removing columns whose sum is 0 (= rectangle with size 0)
    c = each_rect_count[:, each_rect_count.sum(axis=0)>0 ]
    ## compute
    likelihood  = scipy.special.gammaln(alpha * max_elems) * num_blocks
    likelihood -= scipy.special.gammaln(alpha) * max_elems * num_blocks
    likelihood += scipy.special.gammaln(c + alpha).sum()
    likelihood -= scipy.special.gammaln(c.sum(axis=0) + alpha * max_elems).sum()

    return likelihood

def MH_update_input_order(X, rect_locs, row_locs, col_locs, alpha):
    """
    Metropolis-Hastings process for updating the orders of rows and columns in input data

    Parameters
    ----------
    X: numpy.ndarray
        Input data
    rect_locs: numpy.ndarray
        Current rectangular partitioning
    row_locs: numpy.ndarray
        A vector representing relative sorted locations of rows in input data
    col_locs: numpy.ndarray
        A vector representing relative sorted locations of columns in input data
    alpha: float
        A parameter for Dirichlet distribution

    Returns
    ----------
    row_locs: numpy.ndarray
        updated "row_locs"
    col_locs: numpy.ndarray
        updated "col_locs"
    """

    # Initialization
    each_rect_count = count_elems_in_rects(X, rect_locs, row_locs, col_locs)
    prev_likelihood = calc_Dirichlet_likelihood(each_rect_count, alpha)
    num_rows = row_locs.shape[0]
    num_cols = col_locs.shape[0]

    # row-wise Metropolis-Hastings
    l = list(range(num_rows)); random.shuffle(l)
    for ii in l:
        # Computing a likelihood for a rectangle assignment just after modifying one rectangle
        ## backups
        bak = row_locs[ii]
        bak_each_rect_count = each_rect_count.copy()
        ## before updating a row position
        each_rect_count -= count_elems_in_rects(X, rect_locs, row_locs, col_locs, [ii], None)
        ## after updating a row position
        row_locs[ii] = random.uniform()
        each_rect_count += count_elems_in_rects(X, rect_locs, row_locs, col_locs, [ii], None)
        ## computing
        new_likelihood = calc_Dirichlet_likelihood(each_rect_count, alpha)
        # accept or reject
        accept_reject = random.uniform()
        if math.log(accept_reject + np.spacing(accept_reject)) > new_likelihood - prev_likelihood:
            row_locs[ii] = bak
            each_rect_count = bak_each_rect_count
        else:
            prev_likelihood = new_likelihood

    # column-wise Metropolis-Hastings
    l = list(range(num_cols)); random.shuffle(l)
    for ii in l:
        ## backups
        bak = col_locs[ii]
        bak_each_rect_count = each_rect_count.copy()
        ## before updating a row position
        each_rect_count -= count_elems_in_rects(X, rect_locs, row_locs, col_locs, None, [ii])
        ## after updating a row position
        col_locs[ii] = random.uniform()
        each_rect_count += count_elems_in_rects(X, rect_locs, row_locs, col_locs, None, [ii])
        ## computing
        new_likelihood = calc_Dirichlet_likelihood(each_rect_count, alpha)
        # accept or reject: if True, then reject.
        accept_reject = random.uniform()
        if math.log(accept_reject + np.spacing(accept_reject)) > new_likelihood - prev_likelihood:
            col_locs[ii] = bak
            each_rect_count = bak_each_rect_count
        else:
            prev_likelihood = new_likelihood

    return row_locs, col_locs

def MH_update_num_rect(X, bp, beta_rand, rect_locs, uniform_rand, row_locs, col_locs, alpha, enc):
    """
    Metropolis-Hastings process for updating the number of rectangles

    Parameters
    ----------
    X: numpy.ndarray
        Input data
    bp: numpy.ndarray
        A Baxter permutation
    beta_rand: numpy.ndarray
        A random vector drawn from a Beta distribution
    rect_locs: numpy.ndarray
        A matrix representing a rectangular partitioning
    uniform_rand: numpy.ndarray
        A random vector
    row_locs: numpy.ndarray
        A vector representing row locations of rectangles in input data
    col_locs: numpy.ndarray
        A vector representing column locations of rectangles in input data
    alpha: float
        A parameter for Dirichlet distribution
    enc: float
        A desirable number of blocks

    Returns
    ----------
    bp: numpy.ndarray
        updated "bp"
    beta_rand: numpy.ndarray
        updated "beta_rand"
    rect_locs: numpy.ndarray
        updated "rect_locs"
    uniform_rand: numpy.ndarray
        updated "uniform_rand"
    """

    # Initialization
    num_blocks = uniform_rand.shape[0]
    orig_bp = bp.copy()
    orig_beta_rand = beta_rand.copy()
    orig_rect_locs = rect_locs.copy()
    orig_uniform_rand = uniform_rand.copy()

    # Computing a likelihood
    def calc_likelihood(locs, nblocks):
        each_rect_count = count_elems_in_rects(X, locs, row_locs, col_locs)
        likelihood = calc_Dirichlet_likelihood(each_rect_count, alpha)
        likelihood -= ( nblocks - enc * math.log(nblocks) )
        return likelihood

    # the current likelihood
    cur_likelihood  = calc_likelihood(rect_locs, num_blocks)

    # adding or deleting a rectangle
    if random.uniform() > 0.5:
        print('adding a rectangle')
        # adding
        ## lengthening a Baxter permutation
        bp, uniform_rand = evolv_bp(bp, uniform_rand)
        ## converting Baxter permutation to floorplan partitioning
        label_matrix = bp2fp(bp)
        ## creating rectangular partitioning from a given floorplan partitioning
        beta_rand = np.append(beta_rand, random.beta(1, 1))
        rect_locs = evolv_rp(bp, label_matrix, beta_rand[-1], rect_locs)
    elif num_blocks>1:
        print('deleting a rectangle')
        # deleting
        ## converting Baxter permutation to floorplan partitioning
        label_matrix = bp2fp(bp)
        ## deleting a rectangle
        rect_locs = delete_rp(bp, label_matrix, rect_locs)
        ## shortening a Baxter permutation
        bp, uniform_rand = delete_max(bp, uniform_rand)
        beta_rand = np.delete(beta_rand, -1)
    num_blocks = uniform_rand.shape[0]

    # a new likelihood
    new_likelihood  = calc_likelihood(rect_locs, num_blocks)

    # accept or reject: if true, reject
    accept_reject = random.uniform()
    if math.log(accept_reject + np.spacing(accept_reject)) > new_likelihood - cur_likelihood:
        bp = orig_bp
        beta_rand = orig_beta_rand
        rect_locs = orig_rect_locs
        uniform_rand = orig_uniform_rand

    return bp, beta_rand, rect_locs, uniform_rand

def reconstruct_rp(beta_rand, uniform_rand):
    """
    Reconstructing a rectangular partitioning from a given pair of random vectors.

    Parameters
    ----------
    beta_rand: numpy.ndarray
        A random vector drawn from a Beta distribution, controlling rectangle sizes.
    uniform_rand: numpy.ndarray
        A random vector, implicitly representing a Baxter permutation.

    Returns
    ----------
    bp: numpy.ndarray
        A Baxter permutation
    rect_locs: numpy.ndarray
        A matrix representing a rectangular partitioning
    """

    # initialization
    bp = np.array([0])
    rect_locs = np.array([[0, 1, 0, 1]])
    num_blocks = uniform_rand.shape[0]

    # reconstructing rectanglar partitioning step-by-step
    for ii in range(1, num_blocks):
        ## adding an index to a Baxter permutation
        bp, _ = evolv_bp(bp, uniform_rand[:ii], uniform_rand[ii])
        ## converting the Baxter permutation to floorplan partitioning
        label_mat = bp2fp(bp)
        ## creating rectangular partitioning from a given floorplan partitioning
        rect_locs = evolv_rp(bp, label_mat, beta_rand[ii-1], rect_locs)

    return bp, rect_locs

def MH_update_rp(beta_rand, uniform_rand, X, row_locs, col_locs, alpha):
    """
    Metropolis-Hastings updates of rectangular partitioning
    without changing the number of rectangles

    Parameters
    ----------
    beta_rand: numpy.ndarray
        A random vector drawn from a Beta distribution, controlling rectangle sizes
    uniform_rand: numpy.ndarray
        A random vector, implicitly representing a Baxter permutation
    X: numpy.ndarray
        Input data
    row_locs: numpy.ndarray
        A vector representing row locations of rectangles in input data
    col_locs: numpy.ndarray
        A vector representing column locations of rectangles in input data
    alpha: float
        A Dirichlet parameter

    Returns
    ----------
    bp: numpy.ndarray
        A Baxter permutation
    beta_random: numpy.ndarray
        ppdated "beta_random"
    rect_locs: numpy.ndarray
        A matrix representing a rectangular partitioning
    uniform_rand: numpy.ndarray
        updated "uniform_rand"
    """

    # reconstructing a rectangular partitioning from a pair of random variables
    bp, rect_locs = reconstruct_rp(beta_rand, uniform_rand)

    # Computing a likelihood
    def calc_likelihood(locs):
        each_rect_count = count_elems_in_rects(X, locs, row_locs, col_locs)
        likelihood = calc_Dirichlet_likelihood(each_rect_count, alpha)
        return likelihood

    # initial likelihood
    prev_likelihood = calc_likelihood(rect_locs)

    # updating a Baxter permutation
    num_blocks = uniform_rand.shape[0]
    l = list(range(num_blocks)); random.shuffle(l)
    for ii in l:
        # backup
        orig_uniform_rand = uniform_rand[ii]
        # updating
        uniform_rand[ii] = random.uniform()
        new_bp, new_rect_locs = reconstruct_rp(beta_rand, uniform_rand)
        # skipping the following if the Baxter permutation does not change at all
        if (new_bp==bp).all(): continue
        # computing a new likelihood
        new_likelihood = calc_likelihood(new_rect_locs)
        # accept or reject: if true, reject
        accept_reject = random.uniform()
        if math.log(accept_reject + np.spacing(accept_reject)) > new_likelihood - prev_likelihood:
            uniform_rand[ii] = orig_uniform_rand
        else:
            prev_likelihood = new_likelihood
            bp = new_bp
            rect_locs = new_rect_locs

    # updating rectangle sizes
    l = list(range(num_blocks-1)); random.shuffle(l)
    for ii in l:
        # backup
        orig_beta_rand = beta_rand[ii]
        # update
        beta_rand[ii] = random.beta(1, 1)
        _, new_rect_locs = reconstruct_rp(beta_rand, uniform_rand)
        # check
        if (new_rect_locs==rect_locs).all():
            print('skipped b'); continue
        # computing a new likelihood
        new_likelihood = calc_likelihood(new_rect_locs)
        # accept or reject: if true, reject
        accept_reject = random.uniform()
        if math.log(accept_reject + np.spacing(accept_reject)) > new_likelihood - prev_likelihood:
            beta_rand[ii] = orig_beta_rand
        else:
            prev_likelihood = new_likelihood
            rect_locs = new_rect_locs

    return bp, rect_locs, uniform_rand, beta_rand

def calc_test_perplexity(X, rect_locs, row_locs, col_locs, alpha, missing_indices, missing_ratio):
    """
    Computing test perplexity

    Parameters
    ----------
    X: numpy.ndarray
        Input data
    rect_locs: numpy.ndarray
        A matrix representing a rectangular partitioning
    row_locs: numpy.ndarray
        A vector representing row locations of rectangles in input data
    col_locs: numpy.ndarray
        A vector representing column locations of rectangles in input data
    alpha: float
        A Dirichlet parameter
    missing_indices:
        A random matrix with the same size as "X",
        implicitly representing the positions of missing elements used for test data.
    missing_ratio: float
        Ratio of missin elements in the input data

    Returns
    ----------
    perp: float
        Test perplexity
    """

    # Reprecating training data
    X_train = X.copy()
    X_train[ missing_indices < missing_ratio ] = -1

    # Computing Dirichlet pamameters for all the elements
    ### counting elements in blocks
    each_rect_count = count_elems_in_rects(X_train, rect_locs, row_locs, col_locs)
    ## shapes
    max_elems = each_rect_count.shape[0]
    num_blocks = each_rect_count.shape[1]
    ## compute
    each_rect_count_sum = each_rect_count.sum(axis=0)
    keep_blocks = ( each_rect_count_sum>0 )
    norm = np.dot( np.ones((max_elems, 1)), each_rect_count_sum[np.newaxis, keep_blocks] )
    dir_param = (each_rect_count[:, keep_blocks] + alpha) / norm

    # Reprecating test data
    X_test = X.copy()
    X_test[ missing_indices >= missing_ratio ] = -1

    # computing
    each_rect_count = count_elems_in_rects(X_test, rect_locs, row_locs, col_locs)
    ## removing columns whose sum is 0 (= rectangle with size 0)
    cross_entropy = np.sum( each_rect_count[:, keep_blocks] * np.log(dir_param) )

    ## locations of missing elements in the test data
    train_mi = np.zeros(missing_indices.shape, dtype='bool')
    train_mi[ missing_indices < missing_ratio ] = True
    ## locations of originally missing elements
    orig_mi = np.zeros(X.shape, dtype='bool')
    orig_mi[ X<0 ] = True
    ## locations of test data
    test_pos = np.logical_and( train_mi, np.logical_not(orig_mi) )
    ## computing test perplexity
    perp = math.exp( (-1) * cross_entropy / test_pos.sum() )

    return perp

def bpp_mcmc(X, opt=None, verbose=True):
    """
    MCMC sampler for Baxter Permutation process

    Parameters
    ----------
    X: numpy.ndarray
        An input matrix (Currently, we accept matrices, namely 2-dimensional arrays.),
        expected to have dtype=int8 (8-bit signed int).
    opt: utils.Options
        Options, see utils.py for the detail.
    verbose: bool
        Verbose

    Returns
    -------
    rect_locs: numpy.ndarray
        A matrix representing a rectangular partitioning
    row_locs: numpy.ndarray
        A vector representing row locations of rectangles in input data
    col_locs: numpy.ndarray
        A vector representing column locations of rectangles in input data
    perps: numpy.ndarray
        A sequence of test perplexity values
    """

    # option ititialization
    if opt is None:
        opt = utils.BPPOptions(enc=X.size/20000)
    if opt.enc is None:
        opt = utils.BPPOptions(enc=X.size/20000, n=opt.n, alpha=opt.alpha, maxiter=opt.maxiter,
                               missing_ratio=opt.missing_ratio, rand_seed=opt.random_seed)
    if opt.rand_seed is not None:
        np.random.seed(opt.rand_seed)
    else:
        np.random.seed(13)
    if verbose:
        print('The shape of the input data is {}'.format(X.shape))

    # creting test data
    ## keeing the original input data
    original_X = np.copy(X)
    ## masking the input data for testing
    missing_indices = random.uniform(0, 1, X.shape)
    X[ missing_indices < opt.missing_ratio ] = -1

    # Initialization
    if verbose: print('Initialization...')
    ## drawing a random vector from a beta distribution
    beta_rand = random.beta(np.ones(opt.n-1), np.ones(opt.n-1))
    ## creating a Baxter permutation with length "opt.n"
    bp = np.array([0])
    uniform_rand = random.uniform(0, 1, size=1)
    rect_locs = np.array([[0, 1, 0, 1]])  # shape = (1, 4)
    for ii in range(1, opt.n):
        ## adding a new index to the current Baxte permutation
        bp, uniform_rand = evolv_bp(bp, uniform_rand)
        ## converting the Baxter permutation to a floorplan partitioning
        label_matrix = bp2fp(bp)
        ## creating rectangular partitioning from a given floorplan partitioning
        rect_locs = evolv_rp(bp, label_matrix, beta_rand[ii-1], rect_locs)

    # canvas initialization
    figs = plt.subplots(nrows=1, ncols=2, figsize=(12,6))
    plt.ion()  # interactive mode on
    plt.draw(); plt.pause(0.01)  # showing an empty figure

    #
    row_locs = random.uniform(size=X.shape[0])
    col_locs = random.uniform(size=X.shape[1])
    # main loop
    perps = list()
    for iter in range(opt.maxiter):
        if verbose: print('Iter {}...'.format(iter))
        # Metropolis-Hastings update of input data sorting
        if verbose: print('Updating input data sorting...')
        row_locs, col_locs = MH_update_input_order(X, rect_locs, row_locs, col_locs, opt.alpha)
        # Metropolis-Hastings update of the number of rectangles
        if verbose: print('Updating number of rectangles...')
        bp, beta_rand, rect_locs, uniform_rand = \
            MH_update_num_rect(X, bp, beta_rand, rect_locs, uniform_rand, row_locs, col_locs,
                               opt.alpha, opt.enc)
        # Metropolis-Hastings update of rectangles
        if verbose: print('Updating rectangles...')
        bp, rect_locs, uniform_rand, beta_rand = \
            MH_update_rp(beta_rand, uniform_rand, X, row_locs, col_locs, opt.alpha)
        # compute test perplexity
        if verbose: print('rect_locs:\n', rect_locs)
        if verbose: print('Computing test perplexity...')
        perp = calc_test_perplexity(original_X, rect_locs, row_locs, col_locs,
                                    opt.alpha, missing_indices, opt.missing_ratio)
        perps.append(perp)
        # plot
        utils.plot_state(figs, original_X, rect_locs, row_locs, col_locs, perps)

    return rect_locs, row_locs, col_locs, perps
