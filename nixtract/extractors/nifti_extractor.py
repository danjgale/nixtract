
import os
import warnings
import numpy as np
import pandas as pd
import nibabel as nib
from nilearn.input_data import (NiftiMasker, NiftiSpheresMasker, 
                                NiftiLabelsMasker, NiftiMapsMasker)
from nilearn import image
from nilearn.input_data.nifti_spheres_masker import _apply_mask_and_get_affinity

from .base_extractor import BaseExtractor


def _read_coords(roi_file):
    """Parse and validate coordinates from file"""
    if roi_file.endswith('.tsv'):
        coords = pd.read_table(roi_file)
    elif roi_file.endswith('.csv'):
        coords = pd.read_csv(roi_file)
    else:
        raise ValueError('Coordinate file must be a tab-separated .tsv file '
                         'or a comma-separated .csv file')
    
    # validate columns
    required_columns = ['x', 'y', 'z']
    missing_columns = [i for i in required_columns if i not in coords.columns]
    if len(missing_columns) != 0:
        raise ValueError(f'roi_file is missing column headers {missing_columns}')

    # convert to list of lists for nilearn input
    return coords[required_columns].values.tolist()


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

    if not isinstance(roi_file, str):
        raise ValueError('roi_file must be a file name string')

    if roi_file.endswith('.csv') or roi_file.endswith('.tsv'):
        roi = _read_coords(roi_file)
        n_rois = len(roi)
        print('  {} region(s) detected from coordinates'.format(n_rois))

        if kwargs.get('radius') is None:
            warnings.warn('No radius specified for coordinates; setting '
                          'to nilearn.input_data.NiftiSphereMasker default '
                          'of extracting from a single voxel')
        masker = NiftiSpheresMasker(roi, **kwargs)
    
    elif roi_file.endswith('.nii.gz'):
        # remove args for NiftiSpheresMasker 
        if 'radius' in kwargs:
            kwargs.pop('radius')
        if 'allow_overlap' in kwargs:
                kwargs.pop('allow_overlap')
    
        roi_img = image.load_img(roi_file)
        if len(roi_img.shape) == 4:
            n_rois = roi_img.shape[-1]
            print('  {} region(s) detected from {}'.format(n_rois,
                                                        roi_img.get_filename()))
            masker = NiftiMapsMasker(roi_img, allow_overlap=True,**kwargs)
        else:
            n_rois = len(np.unique(roi_img.get_fdata())) - 1
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
    
    else:
        raise ValueError('Invalid file type for roi_file. Must be one of: '
                         '.nii.gz, .csv, .tsv')
    
    return masker, n_rois


class NiftiExtractor(BaseExtractor):
    def __init__(self, fname, roi_file, labels=None, as_voxels=False, 
                 verbose=False, **kwargs):
        """Extract timeseries from a NIFTI image

        Parameters
        ----------
        fname : str
            Functional data
        roi_file : str
            Nifti file containing numeric labels for each voxel to identify
            each region. Can be an atlas/parcellation with multiple regions,
            or a binary mask for a single region. Or, a .tsv file containing
            central coordinates for each region. 
        labels : str, optional
            Label names for each region in roi_file, given in the exact same
            ascending order. If None, a) numeric labels in roi_file will be 
            used if roi_file is a single or multi-region atlas, b) rows in 
            roi_file are enumerated (1-indexed) if roi_file is a .tsv file, or 
            c) voxels are enumerated (1-indexed) if as_voxels is specified. By 
            default None
        as_voxels : bool, optional
            Extract the individual voxel timeseries from a region. Only 
            possible when roi_file is a binary mask (single region), by 
            default False
        verbose : bool, optional
            Print out extraction timestamp, by default False
        **kwargs 
            Arguments to pass to a Nilearn masker object, which is determined
            by the roi_file
        """
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
            return ['voxel{}'.format(int(i))
                    for i in np.arange(self.timeseries.shape[1]) + 1]
        elif isinstance(self.masker, NiftiLabelsMasker): 
            # get actual numerical labels used in image          
            return ['region{}'.format(int(i)) for i in self.masker.labels_]
        elif isinstance(self.masker, NiftiSpheresMasker):
            return ['region{}'.format(int(i)) 
                    for i in np.arange(len(self.masker.seeds)) + 1]
        elif isinstance(self.masker, NiftiMapsMasker):
            return ['region{}'.format(int(i)) 
                    for i in np.arange(self.masker.maps_img.shape[-1]) + 1]

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
        elif isinstance(self.masker, NiftiMapsMasker):
            return self.masker.maps_img_ 