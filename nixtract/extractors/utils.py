
import numpy as np
import pandas as pd
from nilearn import signal


def _mask(darray, roi, as_vertices=False):
    labels = np.unique(roi)
    if len(labels) > 2 and as_vertices:
        raise ValueError('Using as_vertices=True with more than one region '
                         'in roi file. Vertex-level extraction can only be '
                         'performed with a single-region (binary) roi file.')
    if as_vertices:
        timeseries = darray[:, roi.ravel().astype(bool)]
    else:
        timeseries = np.zeros((darray.shape[0], len(labels)))
        for i, l in enumerate(labels):
            mask = np.where(roi == l, 1, 0).astype(bool)
            timeseries[:, i] = darray[:, mask].mean(axis=1)
    
    return timeseries


def mask_data(darray, roi, regressors=None, as_vertices=False, 
              pre_clean=False, **kwargs):
    x = darray.copy()
    if pre_clean:
        x = signal.clean(x, confounds=regressors, **kwargs)
        return _mask(x, roi, as_vertices)
    else:
        timeseries = _mask(x, roi, as_vertices)
        return signal.clean(timeseries, confounds=regressors, **kwargs)


def label_timeseries(tseries, labels, as_vertices):
    if as_vertices:
        cols = [f'vert{i}' for i in np.arange(tseries.shape[1])]
        return pd.DataFrame(tseries, columns=cols)
    else:
        return pd.DataFrame(tseries, columns=labels)