
import os
import numpy as np
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import seaborn as sns
from pandas.plotting import autocorrelation_plot
from nilearn.connectome import sym_matrix_to_vec


def _kde_color():
    return '#404788FF'


def _diverging_palette():
    return sns.diverging_palette(220, 20, s=90, as_cmap=True)


def _check_confounds_alignment(tseries, confounds):
    """Check if tseries has initial scans discarded"""
    diff = confounds.shape[0] - tseries.shape[0]
    if diff > 0:
        confounds = confounds[diff:]
    return confounds


def _title_plot(measures, ax):

    fname = measures['fname']
    confounds = measures['confounds']


    ax.annotate(fname, xy=(-.2, 0.6), xycoords='axes fraction', va='center', 
                ha='left', weight='bold')
    
    ax.annotate(f'confounds: {confounds}', xy=(-.2, 0.3), 
                xycoords='axes fraction', va='center', ha='left')
    plt.axis('off')
    return ax


def _motion_trace_plot(confounds, ax, param_type='trans'):

    list_ = []
    for i in ['x', 'y', 'z']:
        param = f'{param_type}_{i}'
        g = ax.plot(confounds[param], label=i, lw=1)
        list_.append(g)

    xmax = confounds.shape[0]
    interval = xmax / 4
    xticks =  np.round(np.arange(0, xmax + interval, interval), 0).astype(int)
    ax.set_xticks(xticks)
    ax.margins(x=0)

    if param_type == 'trans':
        ax.legend(fontsize='small', frameon=False, labelspacing=0, 
                  ncol=len(list_), bbox_to_anchor=(0, 1.3), 
                  loc='upper left')
        ax.set_ylabel('Translation\n(mm)',  fontsize=8)
    elif param_type == 'rot':
        ax.set_ylabel('Rotation\n(deg)', fontsize=8)

    return ax


def _fd_trace(confounds, ax):

    ax.plot(confounds['framewise_displacement'], c='C3', lw=1)
    ax.axhline(.2, c='k', ls='--')
    ax.margins(x=0)
    ax.set_ylabel('Framewise\n displacement\n(mm)', fontsize=8)

    xmax = confounds.shape[0] - 1
    interval = xmax / 4
    xticks =  np.round(np.arange(0, xmax + interval, interval), 0).astype(int)
    ax.set_xticks(xticks)
    ax.margins(x=0)

    return ax


def _carpet_plot(tseries, fd, ax):

    x = tseries.values.T
    sd = np.mean(x.std(axis=0))
    vmin = x.mean() - (2 * sd)
    vmax = x.mean() + (2 * sd)

    ax.imshow(x, cmap='gray', vmin=vmin, vmax=vmax, aspect='auto')
    ax.set_yticks([])
    ax.set_yticklabels([])

    # set 5 ticks
    xmax = x.shape[1] - 1
    interval = xmax / 4
    xticks =  np.round(np.arange(0, xmax + interval, interval), 0).astype(int)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks + 1)

    ax.set_xlabel('TR')
    ax.set_ylabel('Regions', fontsize=8)

    spikes = np.argwhere(fd > .2).ravel()
    if len(spikes) > 0:
        ypos = ax.get_ylim()[1] - 13
        ax.scatter(x=spikes, y=[ypos] * len(spikes), marker="|", c='C3', s=25, 
                   clip_on=False)

    for d in ["left", "top", "bottom", "right"]:
        plt.gca().spines[d].set_visible(False)
        
    return ax


def _r_plot(mat, ax):
    edges = sym_matrix_to_vec(mat, discard_diagonal=True)
    ax = sns.kdeplot(edges, fill=True, alpha=1, linewidth=0, 
                     color='C7', ax=ax)
    ax.set_xlabel('r')
    ax.set_xlim(-1, 1)
    ax.axvline(0, c='k', ls='--')
    ax.set_ylabel('Density', fontsize=8)

    return ax

def _info_plot(measures, ax):

    mean = measures['mean_r']
    n_sig = measures['sig_edges']
    n_sig_corr = measures['sig_edges_corrected']
    q = measures['q']

    info = (f"Mean r = {mean:.3f}\n"
            f"Prop. p<.05 = {n_sig:.3f}\n"
            f"Prop. q<.05 = {n_sig_corr:.3f}\n"
            f"Modularity (Q) = {q:.3f}")
    ax.annotate(info, xy=(-.5, 0.5), xycoords='axes fraction', va='center', 
                ha='left')
    plt.axis('off')
    return ax


def _tseries_cmat(mat, ax):
    
    ax = sns.heatmap(mat, vmin=-1, vmax=1, cmap=_diverging_palette(), 
                     square=True, ax=ax, cbar_kws={"shrink": .5})
    ax.set_axis_off()
    ax.set_title('Connectivity')
    return ax


def plot_tseries(tseries, confounds, measures, mat, out_dir):
    
    confounds = _check_confounds_alignment(tseries, confounds)

    fig = plt.figure(figsize=(10, 6))
    motion_gs = plt.GridSpec(5, 9, hspace=0.1, wspace=2)
    conn_gs = plt.GridSpec(5, 9, hspace=-.3, wspace=1)

    gs1 = fig.add_subplot(motion_gs[0, :5])
    _title_plot(measures, gs1)

    gs2 = fig.add_subplot(motion_gs[1, :5])
    _motion_trace_plot(confounds, gs2)

    gs3 = fig.add_subplot(motion_gs[2, :5])
    _motion_trace_plot(confounds, gs3, 'rot')

    gs4 = fig.add_subplot(motion_gs[3, :5])
    _fd_trace(confounds, gs4)

    for g in [gs2, gs3, gs4]:
        g.set_xticks([])
        g.set_xticklabels([])

    gs5 = fig.add_subplot(motion_gs[4, :5])
    _carpet_plot(tseries, confounds['framewise_displacement'].values, 
                 gs5)

    gs6 = fig.add_subplot(conn_gs[:4, 5:])
    _tseries_cmat(mat, gs6)

    gs7 = fig.add_subplot(conn_gs[4:, 5:7])
    _r_plot(mat, gs7)

    gs8 = fig.add_subplot(conn_gs[4:, 7])
    _info_plot(measures, gs8)

    out = measures['fname'].replace('timeseries.tsv', 'plot.png')
    fig.savefig(os.path.join(out_dir, out), dpi=300)



def _corrfunc(x, y, ax=None, **kws):
    """Plot the correlation coefficient in the top left hand corner of a plot."""
    r, p = stats.pearsonr(x, y)
    ax = ax or plt.gca()

    xy = (.6, .8) if r < 0 else (.1, .8)
    # sig = 'C3' if p < .05 else 'k'
    # ptext = 'p < .001' if p < .001 else f'p = {p:.3f}'
    
    ax.annotate(f'r = {r:.3f}', xy=xy, fontsize=8, xycoords=ax.transAxes)


def plot_measures(measures, out_dir):
    
    cols = ['mean_fd', 'n_spikes', 'mean_r', 'sig_edges', 'q']
    g = sns.pairplot(measures, vars=cols, diag_kind='kde', 
                     kind='reg', corner=True, height=1.2, aspect=1.7, 
                     plot_kws={'scatter_kws': {'markersize': 5}, 
                               'line_kws': {'c': 'C3'}})
    g.map_lower(_corrfunc)
    g.fig.suptitle('Sample QC measures')

    out = os.path.join(out_dir, 'sample_measures.png')
    g.savefig(out, dpi=300)


def plot_group_connectivity(mat, mean, q, out_dir):

    fig, ax = plt.subplots(figsize=(4, 4))
    fig_pos = ax.get_position()
    # cbar_ax = fig.add_axes([.905, fig_pos.y0, .02, .3])
    ax = sns.heatmap(mat, vmin=-1, vmax=1, cmap=_diverging_palette(), 
                     square=True, ax=ax, cbar_kws={"shrink": .5})
    ax.set_axis_off()

    info = (f"Avg. connectivity matrix\n"
            f"Mean r = {mean:.3f}\n"
            f"Modularity (Q) = {q:.3f}")
    ax.set_title(info)
    out = os.path.join(out_dir, 'connectivity.png')
    fig.savefig(out, dpi=300)


def plot_qc_fc(qc_fc, out_dir):

    fig, ax = plt.subplots(figsize=(3, 2))
    ax = sns.kdeplot(x=qc_fc, fill=True, ax=ax)
    ax.set(xlabel='r')
    ax.axvline(0, c='k', lw=1)
    sns.despine()
    
    # show mean
    med_abs = np.median(np.abs(qc_fc))
    x_loc = ax.get_xlim()[1] * .3
    yloc = ax.get_ylim()[1] * .9   
    ax.text(x=x_loc , y=yloc, s=f'MA = {med_abs:0.3f}')

    fig.subplots_adjust(top=0.9)
    ax.set_title('QC-FC', fontdict={'size': 12})

    out = os.path.join(out_dir, 'qcfc.png')
    fig.savefig(out, dpi=300)


def plot_dist_dependence(distances, qc_fc, r, out_dir):

    g = sns.JointGrid(x=distances, y=qc_fc, ylim=(-1, 1), height=4)
    g.plot_joint(sns.histplot, cmap='viridis', bins=50, cbar=False)
    g.plot_joint(sns.regplot, scatter=False, ci=None, truncate=False, 
                 line_kws=dict(color='m', linewidth=1))
    g.plot_marginals(sns.kdeplot, color=_kde_color(), alpha=.8, linewidth=1, 
                     fill=True)
    g.ax_joint.set(xlabel='Euclidean distance', ylabel='QC-FC', 
                   yticks=np.arange(-1, 1.5, .5))

    # display correlation value
    x_loc = np.max(distances) * .7  
    g.ax_joint.text(x=x_loc , y=.9, s=f'ρ = {r:.3f}', )

    g.fig.subplots_adjust(top=0.9)
    g.fig.suptitle('Distance dependence', fontsize=10)
    
    out = os.path.join(out_dir, 'distance_dependence.png')
    g.fig.savefig(out, dpi=300)
