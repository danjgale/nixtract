
import glob
import natsort
import numpy as np
from scipy.ndimage.measurements import center_of_mass
import pandas as pd
import nibabel as nib

def check_confounds(confounds):
    """Read and verify that required columns are in confounds

    Required columns:
    * 'trans_x'
    * 'trans_y'
    * 'trans_z'
    * 'rot_x'
    * 'rot_y'
    * 'rot_z'
    * 'framewise_displacement'

    These should all be available in fmriprep outputs. 

    Parameters
    ----------
    confounds : str
        File name of confounds file

    Returns
    -------
    pandas.DataFrame
        Loaded confounds

    Raises
    ------
    ValueError
        Any of the required confounds are missing in column headers
    """
    required = ['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', 
                'framewise_displacement']
    try:
        confs = pd.read_table(confounds, usecols=required)
        return confs
    except ValueError as e:
        missing = eval(e.args[0].split(': ')[1])
        raise ValueError(f"{' '.join(missing)} columns not in {confounds}")


def _atlas_to_coords(atlas):
    """Compute center of mass for a volumetric atlas"""
    data = atlas.get_fdata()
    labels = np.unique(data)[1:] # exclude 0
    coords = center_of_mass(data, data, index=labels)
    return np.array(coords)


def read_dist_file(fname):
    """Read file containing distance information

    Parameters
    ----------
    fname : str
        File name of distance file. Can be 1) a .tsv file containing 
        x, y, and z coordinates for each region. If so, columns 'x', 'y', and
        'z' must exist, and rows should included each region. Other columns 
        are ignored. Or, can be 2) a volumetric NIFTI atlas (.nii.gz or .nii)
        in which each region is defined a numeric label. If so, the centers
        of mass will be computed for each region. In either case, rows or 
        labels must be in the same ascending order of the atlas/parcellation 
        used to generate the timeseries.

    Returns
    -------
    numpy.ndarray
        Array of coordinates (n coordinates, 3)

    Raises
    ------
    ValueError
        fname is not read as a nibabel Nifti1Image
    ValueError
        fname does not have column headers 'x', 'y', and 'z'
    ValueError
        fname does not have either .nii, .nii.gz, or .tsv as a file extension
    """
    if fname.endswith('.nii.gz') or fname.endswith('.nii'):
        atlas = nib.load(fname)
        if not isinstance(atlas, nib.Nifti1Image):
            raise ValueError(f'{fname} is not an instance of Nifti1Image')
        coordinates = _atlas_to_coords(atlas)
    elif fname.endswith('.tsv'):
        try:
            coordinates = pd.read_table(fname, usecols=['x', 'y', 'z']).values
        except ValueError as e:
            missing = eval(e.args[0].split(': ')[1])
            raise ValueError(f"{' '.join(missing)} columns not in {fname}")
    else:
        raise ValueError('Distance file must be a NIFTI image (.nii.gz or '
                         '.nii) or a tab-delimited .tsv file')
    
    return coordinates