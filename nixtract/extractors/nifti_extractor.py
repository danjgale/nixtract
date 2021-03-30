
import os
import warnings
import numpy as np
import pandas as pd
import nibabel as nib
from nilearn.input_data import (NiftiMasker, NiftiSpheresMasker, 
                                NiftiLabelsMasker)
from nilearn import image
from nilearn.input_data.nifti_spheres_masker import _apply_mask_and_get_affinity

from .base_extractor import BaseExtractor


def _read_coords(roi_file):
    """Parse and validate coordinates from file"""

    if not roi_file.endswith('.tsv'):
        raise ValueError('Coordinate file must be a tab-separated .tsv file')

    coords = pd.read_table(roi_file)
    
    # validate columns
    columns = [x for x in coords.columns if x in ['x', 'y', 'z']]
    if (len(columns) != 3) or (len(np.unique(columns)) != 3):
        raise ValueError('Provided coordinates do not have 3 columns with '
                         'names `x`, `y`, and `z`')

    # convert to list of lists for nilearn input
    return coords.values.tolist()


def _get_spheres_from_masker(masker, ref_img):
    """Re-extract spheres from coordinates to make niimg. 
    
    Note that this will take a while, as it uses the exact same function that
    nilearn calls to extract data for NiftiSpheresMasker
    """
    ref_img = nib.Nifti1Image(ref_img.get_fdata()[:, :, :, [0]], 
                              ref_img.affine)

    X, A = _apply_mask_and_get_affinity(masker.seeds, ref_img, masker.radius, 
                                        masker.allow_overlap)
    # label sphere masks
    spheres = A.toarray()
    spheres *= np.arange(1, len(masker.seeds) + 1)[:, np.newaxis]

    # combine masks, taking the maximum if overlap occurs
    arr = np.zeros(spheres.shape[1])
    for i in np.arange(spheres.shape[0]):
        arr = np.maximum(arr, spheres[i, :])
    arr = arr.reshape(ref_img.shape[:-1])
    spheres_img = nib.Nifti1Image(arr, ref_img.affine)
    
    if masker.mask_img is not None:
        mask_img_ = image.resample_to_img(masker.mask_img, spheres_img)
        spheres_img = image.math_img('img1 * img2', img1=spheres_img, 
                               img2=mask_img_)
    return spheres_img


def _set_volume_masker(roi_file, as_voxels=False, **kwargs):
    """Check and see if multiple ROIs exist in atlas file"""

    if isinstance(roi_file, str) and roi_file.endswith('.tsv'):
        roi = _read_coords(roi_file)
        n_rois = len(roi)
        print('  {} region(s) detected from coordinates'.format(n_rois))

        if kwargs.get('radius') is None:
            warnings.warn('No radius specified for coordinates; setting '
                          'to nilearn.input_data.NiftiSphereMasker default '
                          'of extracting from a single voxel')
        masker = NiftiSpheresMasker(roi, **kwargs)
    else:
        # remove args for NiftiSpheresMasker 
        if 'radius' in kwargs:
            kwargs.pop('radius')
        if 'allow_overlap' in kwargs:
            kwargs.pop('allow_overlap')
    
        roi_img = image.load_img(roi_file)
        n_rois = len(np.unique(roi_img.get_data())) - 1
        print('  {} region(s) detected from {}'.format(n_rois,
                                                       roi_img.get_filename()))
        if n_rois > 1:
            masker = NiftiLabelsMasker(roi_img, **kwargs)
        elif n_rois == 1:
            # binary mask for single ROI 
            if as_voxels:
                if 'mask_img' in kwargs:
                    kwargs.pop('mask_img')
                masker = NiftiMasker(roi_img, **kwargs)
            else:
                # more computationally efficient if only wanting the mean
                masker = NiftiLabelsMasker(roi_img, **kwargs)
        else:
            raise ValueError('No ROI detected; check ROI file')
    
    return masker, n_rois


class NiftiExtractor(BaseExtractor):
    def __init__(self, fname, roi_file, labels=None, as_voxels=False, 
                 verbose=False, **kwargs):

        self.fname = fname
        self.img = nib.load(fname)
        self.roi_file = roi_file
        self.labels = labels
        self.as_voxels = as_voxels
        self.verbose = verbose

        # determine masker
        self.masker, self.n_rois = _set_volume_masker(roi_file, as_voxels, 
                                                      **kwargs)
        self.masker_type = self.masker.__class__.__name__
        self.regressor_names = None
        self.regressor_array = None
        
    def _get_default_labels(self):
        """Generate default numerical (1-indexed) labels depending on the 
        masker
        """
        self.check_extracted()
        
        if isinstance(self.masker, NiftiMasker):
            return ['voxel {}'.format(int(i))
                    for i in np.arange(self.data.shape[1]) + 1]
        elif isinstance(self.masker, NiftiLabelsMasker): 
            # get actual numerical labels used in image          
            return ['roi {}'.format(int(i)) for i in self.masker.labels_]
        elif isinstance(self.masker, NiftiSpheresMasker):
            return ['roi {}'.format(int(i)) 
                    for i in np.arange(len(self.masker.seeds)) + 1]

    def discard_scans(self, n_scans):
        """Discard first N scans from data and regressors, if available 

        Parameters
        ----------
        n_scans : int
            Number of initial scans to remove
        """
        arr = self.img.get_data()
        arr = arr[:, :, :, n_scans:]
        self.img = nib.Nifti1Image(arr, self.img.affine)

        if self.regressor_array is not None:
            self.regressor_array = self.regressor_array[n_scans:, :]
        
        return self

    def extract(self):
        """Extract timeseries data using the determined nilearn masker"""
        self.show_extract_msg(self.fname)
        timeseries = self.masker.fit_transform(self.img, 
                                               confounds=self.regressor_array)
        self.timeseries = pd.DataFrame(timeseries)
        
        if self.labels is None:
            self.timeseries.columns = self._get_default_labels()
        else:
            self.timeseries.columns = self.labels
        
        return self

    def get_fitted_roi_img(self):
        """Return fitted roi img from nilearn maskers

        Returns
        -------
        nibabel.Nifti1Image
            Image generated and used by underlying nilearn masker class.  
        """
        if isinstance(self.masker, NiftiMasker):
            return self.masker.mask_img_
        elif isinstance(self.masker, NiftiLabelsMasker):
            return self.masker.labels_img
        elif isinstance(self.masker, NiftiSpheresMasker):
            return _get_spheres_from_masker(self.masker, self.img)