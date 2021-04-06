
import os
import warnings
from itertools import repeat
import multiprocessing
from datetime import datetime
import numpy as np
from scipy import stats
from scipy.spatial.distance import cdist
import pandas as pd
from nilearn.connectome import ConnectivityMeasure, sym_matrix_to_vec
import bct

from .utils import check_confounds


def _r_to_p(r, n):
    """Convert a Pearson r value to p values"""
    t = (r * np.sqrt(n - 2)) / np.sqrt(1 - r**2)
    return stats.t.sf(t, df=n - 2) * 2


# def _fdr(pvals):
#     """Benjamini-Hochberg false-discovery rate correction

#     Note: This function from Pingouin, but is adapted and implemented here to 
#     avoid requiring an entire dependency for just one function. This function 
#     is specific to QC analyses and not intended for general use. The original 
#     code can be found at:
#     https://github.com/raphaelvallat/pingouin/blob/91b184780d27f291ff9f894d1c80f74c0c20e146/pingouin/multicomp.py#L12-L117

#     Please use the Pinguoin function if you wish to include FDR corrections in 
#     your own separate analyses.

#     Parameters
#     ----------
#     pvals : array_like
#         Array of p-values of the individual tests.

#     Returns
#     -------
#     pval_corrected : array
#         Corrected p values from FDR correction
#     """
#     # Convert to array and save original shape
#     pvals = np.asarray(pvals)
#     shape_init = pvals.shape
#     pvals = pvals.ravel()
#     num_nan = np.isnan(pvals).sum()

#     # Sort the (flattened) p-values
#     pvals_sortind = np.argsort(pvals)
#     pvals_sorted = pvals[pvals_sortind]
#     sortrevind = pvals_sortind.argsort()
#     ntests = pvals.size - num_nan

#     # Empirical CDF factor
#     ecdffactor = np.arange(1, ntests + 1) / float(ntests)

#     # Now we adjust the p-values
#     pvals_corr = np.diag(pvals_sorted / ecdffactor[..., None])
#     pvals_corr = np.minimum.accumulate(pvals_corr[::-1])[::-1]
#     pvals_corr = np.clip(pvals_corr, None, 1)

#     # And revert to the original shape and order
#     pvals_corr = np.append(pvals_corr, np.full(num_nan, np.nan))
#     pvals_corrected = pvals_corr[sortrevind].reshape(shape_init)

#     return pvals_corrected


def count_sig_edges(x, n):
    """Count proportion of significant edges in functional connectivity matrix
    
    Proportion of p values before and after FDR correction are both computed 

    Parameters
    ----------
    x : numpy.ndarray
        Dense connectivity matrix of Pearson r values
    n : int
        Number of timepoints in timeseries

    Returns
    -------
    float, float
        Proportion of significant edges, uncorrected and FDR-corrected, 
        respectively
    """
    edges = sym_matrix_to_vec(x, discard_diagonal=True)
    pvals = np.array([_r_to_p(r, n) for r in edges])
    prop_sig = len(pvals[pvals < .05]) / len(pvals)

    # fdr_pvals = _fdr(pvals)
    # corrected_prop_sig = len(fdr_pvals[fdr_pvals < .05]) / len(fdr_pvals)

    return prop_sig#, corrected_prop_sig


def network_modularity(x, n_iters=100):
    """Mean network modularity, Q, estimated using repeated iterations of
    the Louvain algorithm (modified for positive and negative weights)

    Parameters
    ----------
    x : numpy.ndarray
        Dense connectivity matrix of Pearson r values
    n_iters : int, optional
        Number of iterations, by default 100

    Returns
    -------
    float
        The mean Q statistic across iterations
    """
    modularities = []
    for i in range(n_iters):
        _, q = bct.modularity_louvain_und_sign(x)
        modularities.append(q)
    return np.mean(modularities)
    

def count_spikes(x, thresh):
    """Count the number of spikes above a threshold for a given confound

    Parameters
    ----------
    x : array-like
        Timeseries of confound of interest (e.g., framewise_displacement)
    thresh : float
        Spike threshold

    Returns
    -------
    int
        Number of timepoints that exceed threshold
    """
    return len(x[x > thresh])


def qc_fc(x, fd):
    """Correlate Pearson r and framewise displacement for each edge in 
    connectivity matrix

    Parameters
    ----------
    x : list, (n timeseries)
        List of dense connectivity matrices using Pearson r values
    fd : array-like, (n timeseries)
        The mean framewise displacement for associated with each connectivity 
        matrix in `x`

    Returns
    -------
    list
        The Spearman r correlations for each edge
    """
    edges = np.array([sym_matrix_to_vec(i, True) for i in x])
    return [stats.spearmanr(edges[:, i], fd)[0] for i in np.arange(edges.shape[1])]


def compute_tseries_measures(tseries, confounds):
    """Compute quality measures for a single timeseries dataframe

    Measures: 
    * Mean framewise displacement
    * Number of framewise displacement spikes
    * Unthresholded connectivity matrix 
    * Mean connectivity
    * Number of significant edges in connectivity matrix (uncorrected and FDR-
      corrected)
    * Modularity, Q

    Parameters
    ----------
    tseries : pandas.DataFrame, (n timepoints, n signals)
        Timeseries data where each column is a separate signal (i.e. region) 
        and rows are timepoints
    confounds : pandas.DataFrame (n timepoints, n confounds)
        Confounds data where each column is a separate confound and rows are
        timepoints 

    Returns
    -------
    dict, numpy.ndarray
        Quality metrics and connectivity matrix
    """
    n_samples = tseries.shape[0]
    n_spikes = count_spikes(confounds['framewise_displacement'], thresh=.2)
    mean_fd = confounds['framewise_displacement'].mean()
    
    # connectivity
    cm = ConnectivityMeasure(kind='correlation')
    mat = cm.fit_transform([tseries.values])[0]

    mean_r = np.mean(sym_matrix_to_vec(mat, discard_diagonal=True))
    count = count_sig_edges(mat, n_samples)
    modularity = network_modularity(mat)
    
    measures = {
        'n': n_samples,
        'mean_fd': mean_fd,
        'n_spikes': n_spikes,
        'mean_r': mean_r,
        'sig_edges': count, 
        # 'significant_edges_corrected': corrected_count,
        'q': modularity
    }
    return measures, mat


def analyze_tseries(fname, confounds, plot=True, out_dir=None, verbose=False):
    """Perform quality analysis on a single timeseries file

    Parameters
    ----------
    fname : str
        Path to timeseries file
    confounds : str
        Path to confounds file, must contain columns headers 'trans_x', 
        'trans_y' 'trans_z', 'rot_x', 'rot_y', 'rot_z', and
        'framewise_displacement'
    plot : bool, optional
        Generate a summary plot, by default True
    out_dir : str, optional
        Path to save summary plot if specified. No plot will be saved if None. 
        By default None

    Returns
    -------
    dict, numpy.ndarray
        Quality metrics and connectivity matrix
    """
    if verbose:
        t = datetime.now().strftime("%H:%M:%S")
        print(f'[{t}] Analyzing {os.path.basename(fname)}')
    tseries = pd.read_table(fname)
    confounds_df = check_confounds(confounds)
    measures, mat = compute_tseries_measures(tseries, confounds_df)
    measures['fname'] = os.path.basename(fname)
    measures['confounds'] = os.path.basename(confounds)
    # if plot and save_dir:
    #     plot_scan(tseries, confounds, measures, mat, out_dir)
    return measures, mat


def compute_dataset_measures(fc_matrices, measures, out_dir, coords=None):
    """Compute group-level quality measures

    Measures: 
    * Unthresholded group average connectivity matrix
    * Mean connectivity of group average connectivity matrix
    * Number of significant edges in the group average connectivity matrix 
      (uncorrected and FDR-corrected), 
    * Modularity, Q
    * QC-FC
    * Distance dependance QC-FC (if coordinates are provided)

    Parameters
    ----------
    fc_matrices : list (n timeseries)
        List of dense connectivity matrices using Pearson r values
    measures : pandas.DataFrame (n timeseries, n measures)
        Dataframe with quality measures for each timeseries (i.e. the output
        from `analyze_tseries`)
    out_dir : str
        Output directory to save plots
    coords : numpy.array, optional (n regions, 3)
        Coordinate array to compute distances between each region in timeseries,
        which is required for distance dependance QC-FC. Columns should be
        X, Y, and Z, respectively. By default None
    """
    # average fc matrix and measures (sig. measures not performed as r 
    # values in matrix don't reflect actual corr stats)
    group_mean_fc = np.mean(fc_matrices, axis=0)
    mean_fc = np.mean(sym_matrix_to_vec(group_mean_fc, discard_diagonal=True))
    modularity_abs = network_modularity(group_mean_fc)
    modularity_prop = network_modularity(group_mean_fc)
    # group_connectivity_plot(group_mean_fc, mean_fc, modularity_abs, 
    #                         modularity_prop, out_dir)

    # QC-FC
    qcfc_data = qc_fc(fc_matrices, measures['mean_fd'])
    median_abs_qcfc = np.median(np.abs(qcfc_data))
    # qcfc_plot() 
    
    if coords is not None:
        distances = sym_matrix_to_vec(cdist(coords, coords), 
                                      discard_diagonal=True)
        dist_dependence = stats.spearmanr(qcfc_data, distances)
        # dist_dependence_plot()
    else:
        print('No atlas provided, skipping distance dependence QC-FC')


def quality_analysis(timeseries, confounds, coords, out_dir, group_only=False,
                     n_jobs=1, verbose=False):
    """Perform full quality analysis on a timeseries dataset

    Measures for each timeseries: 
    * Mean framewise displacement
    * Number of framewise displacement spikes
    * Unthresholded connectivity matrix 
    * Mean connectivity
    * Number of significant edges in connectivity matrix (uncorrected and FDR-
      corrected)
    * Modularity, Q

    Individual timeseries measures are saved in `measures.tsv` in `out_dir`. 

    Group-level measures across timeseries: 
    * Unthresholded group average connectivity matrix
    * Mean connectivity of group average connectivity matrix
    * Number of significant edges in the group average connectivity matrix 
      (uncorrected and FDR-corrected), 
    * Modularity, Q
    * QC-FC
    * Distance dependance QC-FC (if coordinates are provided)

    Note that group-level measures are skipped if only one timeseries file is 
    provided. As well, group-level measures may not be stable with small 
    sample sizes.

    Parameters
    ----------
    timeseries : list of str
        File names for each timeseries to be analyzed
    confounds : list of str
        File names for confounds files associated with each timeseries. Should
        be in the exact same order.
    coords : numpy.array, optional (n regions, 3)
        Coordinate array to compute distances between each region in timeseries,
        which is required for distance dependance QC-FC. Columns should be
        X, Y, and Z, respectively. By default None
    out_dir : str
        Output directory. Plots are saved in out_dir/plots
    group_only : bool, optional
        Skip plots for individual timeseries and only generate group-level 
        plots. Recommended is time or file limits are a concern. By default 
        False
    n_jobs : int, optional
        The number of CPUs to use if parallelization is desired, by default 1
        (serial processing)

    Raises
    ------
    ValueError
        Timeseries and confounds are different lengths
    ValueError
        Timeseries and confounds are empty lists
    """
    if len(timeseries) == len(confounds):
        n_ts = len(timeseries)
    else:
        raise ValueError(f'The number of timeseries files ({len(timeseries)}) '
                         'does not equal the number of confounds files '
                         f'({len(confounds)})')
    if n_ts == 0:
        raise ValueError('No timeseries and confounds files provided')
    if n_ts < 10:
        warnings.warn('Fewer than 10 timeseries files detected! Group-level '
                      'measures (e.g., QC-FC) may not be stable')

    os.makedirs(out_dir, exist_ok=True)
    plot_scans = True if not group_only else False
    plot_dir = os.path.join(out_dir, 'plots')
    os.makedirs(plot_dir, exist_ok=True)
    
    # scan-level measures and plots
    if n_jobs == 1:
        ts_measures, fc_matrices = [], []
        for ts, confs in zip(timeseries, confounds):
            data, res, mat = analyze_tseries(ts, confs, plot=plot_scans, 
                                             out_dir=plot_dir, verbose=verbose)
            ts_measures.append(res)
            fc_matrices.append(mat)

    else:
        args = zip(timeseries, confounds, repeat(plot_scans), repeat(plot_dir),
                   repeat(verbose))
        with multiprocessing.Pool(processes=n_jobs) as pool:
            res = pool.starmap(analyze_tseries, args)
        ts_measures = [i[0] for i in res]
        fc_matrices = [i[1] for i in res]

    measures = pd.DataFrame(ts_measures)
    measures = measures[['fname', 'confounds', 'n', 'mean_fd', 'n_spikes', 
                         'mean_r', 'sig_edges', 'q']]
    measures.to_csv(os.path.join(out_dir, 'measures.tsv'), sep='\t')

    # group-level measures
    if n_ts > 1:
        if verbose:
            t = datetime.now().strftime("%H:%M:%S")
            print(f'[{t}] Computing group-level measures')
        compute_dataset_measures(fc_matrices, measures, out_dir, 
                                 coords=coords)
    else:
        warnings.warn('Only one timeseries file provided; skipping group '
                      'measures')
