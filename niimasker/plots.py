"""Generate figures for visual report"""

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from nilearn.plotting import plot_roi, plot_matrix
from nilearn.image import mean_img
from nilearn.connectome import ConnectivityMeasure
from niimasker.report import make_report


def plot_timeseries(data, cmap):
    """Plot timeseries traces for each extracted ROI.

    Parameters
    ----------
    data : pandas.core.DataFrame
        Timeseries data extracted from niimasker.py
    cmap : matplotlib.colors.LinearSegmentedColormap
        Colormap to use.
    Returns
    -------
    matplotlib.pyplot.figure
        Timeseries plot
    """
    n_rois = data.shape[1]
    fig, axes = plt.subplots(n_rois, 1, sharex=True,
                             figsize=(15, int(n_rois / 5)))

    cmap_vals = np.linspace(0, 1, num=n_rois)

    for i in np.arange(n_rois):

        ax = axes[i]
        y = data.iloc[:, i]
        x = y.index.values

        # draw plot
        ax.plot(x, y, c=cmap(cmap_vals[i]))
        ax.set_ylabel(data.columns[i], rotation='horizontal',
                      position=(-.1, -.1), ha='right')

        # remove axes and ticks
        plt.setp(ax.spines.values(), visible=False)
        ax.tick_params(left=False, labelleft=False)
        ax.xaxis.set_visible(False)

    fig.tight_layout()
    return fig


def plot_carpet(data):

    plot_data = data.transpose().values
    vlim = np.max(np.abs(plot_data))
    fig, axes = plt.subplots(2, 1, sharex=True, figsize=(15, 8),
                           gridspec_kw={'height_ratios': [.2, 1]})
    axes[0].plot(np.arange(plot_data.shape[1]), np.mean(plot_data, axis=0),
                 c='k')
    axes[0].set_ylabel('Mean BOLD')
    im = axes[1].imshow(plot_data, cmap='coolwarm', aspect='auto', vmin=-vlim,
                   vmax=vlim)
    cbar = axes[1].figure.colorbar(im, ax=axes[1], orientation='horizontal',
                                   fraction=.05)
    axes[1].set_ylabel('Voxels')
    axes[1].set_xlabel('Volumes')
    fig.tight_layout()
    return fig


def plot_mask(mask_img, func_img, cmap):
    """Overlay mask/atlas on mean functional image.

    Parameters
    ----------
    atlas_img : str
        File name of atlas/mask image
    func_img : str
        File name of 4D functional image that was used in extraction.
    cmap : matplotlib.colors.LinearSegmentedColormap
        Colormap to use.

    Returns
    -------
    matplotlib.pyplot.figure
        Atlas/mask plot
    """
    # compute mean of functional image
    bg_img = mean_img(func_img)

    n_cuts = 7
    fig, axes = plt.subplots(3, 1, figsize=(15, 6))

    g = plot_roi(mask_img, bg_img=bg_img, display_mode='z', axes=axes[0],
                 alpha=.66, cut_coords=np.linspace(-50, 60, num=n_cuts),
                 cmap=cmap, black_bg=True, annotate=False)
    g.annotate(size=8)
    g = plot_roi(mask_img,  bg_img=bg_img, display_mode='x', axes=axes[1],
                 alpha=.66, cut_coords=np.linspace(-60, 60, num=n_cuts),
                 cmap=cmap, black_bg=True, annotate=False)
    g.annotate(size=8)
    g = plot_roi(mask_img, bg_img=bg_img, display_mode='y', axes=axes[2],
                 alpha=.66, cut_coords=np.linspace(-90, 60, num=n_cuts),
                 cmap=cmap, black_bg=True, annotate=False)
    g.annotate(size=8)
    return fig


def plot_connectome(data, tick_cmap, labels=None):

    cm = ConnectivityMeasure(kind='correlation')
    mat = cm.fit_transform([data])[0]

    labels = [u"\u25A0"] * data.shape[1]

    fig, ax = plt.subplots(figsize=(15, 15))
    plot_matrix(mat, labels=labels, tri='lower', figure=fig, vmin=-1, vmax=1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    cmap_vals = np.linspace(0, 1, num=len(labels))
    for i, lab in enumerate(labels):
        ax.get_xticklabels()[i].set_color(tick_cmap(cmap_vals[i]))
        ax.get_yticklabels()[i].set_color(tick_cmap(cmap_vals[i]))
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontdict={'fontsize':10, 'verticalalignment': 'center', 'horizontalalignment': 'center'})
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontdict={'fontsize':10, 'verticalalignment': 'center', 'horizontalalignment': 'center'})
    return fig


def plot_qcfc(motion_metric):
    pass


def make_figures(functional_images, timeseries_dir, mask_img, as_carpet=False,
                 connectivity_metrics=True, motion_metric=None):

    figure_dir = os.path.join(timeseries_dir, 'niimasker_data/figures/')
    os.makedirs(figure_dir, exist_ok=True)
    report_dir = os.path.join(timeseries_dir, 'reports')
    os.makedirs(report_dir, exist_ok=True)

    for func in functional_images:

        func_img_name = os.path.basename(func).split('.')[0]

        timeseries_file = os.path.join(timeseries_dir,
                                       '{}_timeseries.tsv'.format(func_img_name))
        timeseries_data = pd.read_csv(timeseries_file, sep=r'\t', engine='python')

        n_rois = timeseries_data.shape[1]
        if n_rois > 1:
            roi_cmap = matplotlib.cm.get_cmap('binary')
        else:
            roi_cmap = matplotlib.cm.get_cmap('nipy_spectral')

        # plot and save timeseries
        if as_carpet:
            fig = plot_carpet(timeseries_data)
            bbox_inches = 'tight'
        else:
            fig = plot_timeseries(timeseries_data, roi_cmap)
            bbox_inches = None
        timeseries_fig = os.path.join(figure_dir,
                                      '{}_timeseries_plot.png'.format(func_img_name))
        fig.savefig(timeseries_fig, bbox_inches=bbox_inches)
        # plot and save mask overlay
        fig = plot_mask(mask_img, func, roi_cmap)
        overlay_fig = os.path.join(figure_dir,
                                 '{}_atlas_plot.png'.format(func_img_name))
        fig.savefig(overlay_fig, bbox_inches='tight')

        # place-holder for connectivity plots
        if connectivity_metrics:
            fig = plot_connectome(timeseries_data.values, roi_cmap,
                                  timeseries_data.columns)
            connectome_fig = os.path.join(figure_dir, '{}_connectome_plot.png'.format(func_img_name))
            fig.savefig(connectome_fig, bbox_inches='tight')
            # plot_qcfc(motion_metric)
            qcfc_fig = os.path.join(figure_dir, '{}_qcfc_plot.png'.format(func_img_name))
            # fig.savefig(qcfc_fig)
        else:
            connectome_fig = None
            qcfc_fig = None

        # generate report
        make_report(func, timeseries_dir, overlay_fig, timeseries_fig,
                    connectome_fig, qcfc_fig)
