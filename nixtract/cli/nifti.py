"""Functions for command line interface
"""
import argparse
import os
import shutil
import pandas as pd
from nilearn.datasets import (fetch_atlas_destrieux_2009, fetch_atlas_yeo_2011,
                              fetch_atlas_aal, fetch_atlas_basc_multiscale_2015,
                              fetch_atlas_talairach, fetch_atlas_schaefer_2018)

from nixtract.cli.base import (base_cli, handle_base_args, replace_file_ext,
                               make_param_file, check_glob, empty_to_none, 
                               run_extraction)
from nixtract.extractors import NiftiExtractor

def _cli_parser():
    """Reads NIFTI CLI arguments and returns input specifications combined with
    those from the general CLI
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_files', nargs='+', type=str,
                        help='One or more input NIFTI images (.nii.gz only). '
                             'Can also be a single string with wildcards (*) '
                             'to specify all files matching the file pattern. '
                             'If so, these files are naturally sorted by file '
                             'name prior to extraction')
    parser.add_argument('--roi_file', type=str, 
                        help='Parameter that defines the region(s) of interest. '
                             'This can be 1) a file path to NIFTI image that is '
                             'an atlas of multiple regions, a probabilistic atlas, '
                             'or a binary mask of one region, 2) a nilearn query string '
                             'formatted as `nilearn:<atlas-name>:<atlas-parameters> '
                             '3) a file path to a .tsv file that has x, y, z columns '
                             'that contain coordinates in MNI space. Refer to '
                             'online documentation for more detail and how '
                             'these options map onto the underlying nilearn '
                             'masker classes')
    parser.add_argument('--mask_img', type=str,
                        help='File path of a binary mask a to be used when '
                             '`roi_file` is a) an multi-region atlas or b) a list '
                             'of coordinates. This will restrict extraction to '
                             'only voxels within the mask. If `roi_file` is a '
                             'single region binary mask, this will be ignored')
    parser.add_argument('--labels', nargs='+', type=str, 
                        help='Labels corresponding to the region numbers in '
                             '`roi_file`. Can either be a) a list of strings, b) '
                             'or a .tsv file that contains a `Labels` column. '
                             'Labels must be sorted in ascending order to correctly '
                             'correspond to the atlas indices. The number of '
                             'labels provided must match the number of non-zero '
                             'indices in `roi_file`. Numeric indices are used if '
                             'not provided')
    parser.add_argument('--as_voxels', default=False, action='store_true',
                        help='Extract the timeseries of each voxel instead in a '
                             'a region rather than the mean timeseries. This is '
                             'only available for single region (binary) masks. '
                             'Default: False')
    parser.add_argument('--radius', type=float, metavar='radius', 
                        help='Set the radius of the spheres (in mm) centered on '
                             'the coordinates provided in `roi_file`. Only applicable '
                             'when a coordinate .tsv file is passed to `roi_file`; '
                             'otherwise, this will be ignored. If not set, '
                             'the timeseries of each coordinate is extracted '
                             '(nilearn default)')
    parser.add_argument('--allow_overlap', action='store_true', default=False,
                        help='Permit overlapping spheres when coordinates are '
                             'provided to `roi_file` and `radius` is provided')
    parser.add_argument('--smoothing_fwhm', type=float,
                        help='Smoothing kernel FWHM (in mm) if spatial smoothing '
                             'is desired.')
    parser = base_cli(parser)
    return parser.parse_args()


def get_labelled_atlas(query, data_dir=None, return_labels=True):
    """Parses input query to determine which atlas to fetch and what version
    of the atlas to use (if applicable).

    Parameters
    ----------
    query : str
        Input string in the following format:
        nilearn:{atlas_name}:{atlas_parameters}. The following can be for
        `atlas_name`: 'destrieux', 'yeo', 'aal', 'talairach', and 'schaefer'.
        `atlas_parameters` is not available for the `destrieux` atlas.
    data_dir : str, optional
        Directory in which to save atlas data. By default None, which creates
        a ~/nilearn_data/ directory as per nilearn.
    return_labels : bool, optional
        Whether to return atlas labels. Default is True. Not available for the
        'basc' atlas.

    Returns
    -------
    str, list or None
        The atlas image and the accompanying labels (if provided)

    Raises
    ------
    ValueError
        Raised when the query does is not formatted correctly or if the no
        match found.
    """
    # extract parameters
    params = query.split(':')
    if len(params) == 3:
        _, atlas_name, sub_param = params
    elif len(params) == 2:
        _, atlas_name = params
        sub_param = None
    else:
        raise ValueError('Incorrect atlas query string provided')

    # get atlas
    if atlas_name == 'destrieux':
        atlas = fetch_atlas_destrieux_2009(lateralized=True, data_dir=data_dir)
        img = atlas['maps']
        labels = atlas['labels']
    elif atlas_name == 'yeo':
        atlas = fetch_atlas_yeo_2011(data_dir=data_dir)
        img = atlas[sub_param]
        if '17' in sub_param:
            labels = pd.read_csv(atlas['colors_17'], sep=r'\s+')['NONE'].tolist()
    elif atlas_name == 'aal':
        version = 'SPM12' if sub_param is None else sub_param
        atlas = fetch_atlas_aal(version=version, data_dir=data_dir)
        img = atlas['maps']
        labels = atlas['labels']
    elif atlas_name == 'basc':

        version, scale = sub_param.split('-')
        atlas = fetch_atlas_basc_multiscale_2015(version=version,
                                                 data_dir=data_dir)
        img = atlas['scale{}'.format(scale.zfill(3))]
        labels = None
    elif atlas_name == 'talairach':
        atlas = fetch_atlas_talairach(level_name=sub_param, data_dir=data_dir)
        img = atlas['maps']
        labels = atlas['labels']
    elif atlas_name == 'schaefer':
        n_rois, networks, resolution = sub_param.split('-')
        atlas = fetch_atlas_schaefer_2018(n_rois=int(n_rois),
                                          yeo_networks=int(networks),
                                          resolution_mm=int(resolution),
                                          data_dir=data_dir)
        img = atlas['maps']
        labels = atlas['labels']
    else:
        raise ValueError('No atlas detected. Check query string')

    if not return_labels:
        labels = None
    else:
        labels = labels.astype(str).tolist()

    return img, labels


def _check_nifti_params(params):
    """Ensure that required fields are included and correctly formatted"""
    params = handle_base_args(params)

    if params['input_files'] is None:
        raise ValueError('Missing input files. Check files')
    else:
        params['input_files'] = check_glob(params['input_files'])
        # glob returned nothing
        if not params['input_files']:
            raise ValueError('Missing input files. Check files')

    if not all([i.endswith('.nii.gz') for i in params['input_files']]):
        raise ValueError('input_files must be gzipped NIFTI images')

    if not params['roi_file']:
        raise ValueError('Missing roi_file input.')
    
    if params['roi_file'].startswith('nilearn:'):
        cache = os.path.join(params['output_dir'], 'nixtract_data')
        os.makedirs(cache, exist_ok=True)
        atlas, labels = get_labelled_atlas(params['roi_file'], data_dir=cache,
                                           return_labels=True)
        params['roi_file'] = atlas
        params['labels'] = labels

    params['labels'] = empty_to_none(params['labels'])
    if isinstance(params['labels'], str):
        if params['labels'].endswith('.tsv'):
            df = pd.read_table(params['labels'])
            params['labels'] = df['Label'].tolist()
        else:
            raise ValueError('Labels must be a filename or a list of strings.')

    return params


def extract_nifti(input_file, roi_file, regressor_file, params):
    """Extract timeseries from a NIFTI image

    Parameters
    ----------
    input_files : str
        File path of the functional nifti file
    roi_file : str
        File path of the roi file, where each region is labelled based on the
        numeric values in the file 
    regressor_file : str
        File path of regressor file
    params : dict
        Parameter dictionary for extraction
    """
    extractor = NiftiExtractor(
        fname=input_file,
        roi_file=roi_file, 
        labels=params['labels'],
        as_voxels=params['as_voxels'],
        mask_img=params['mask_img'], 
        radius=params['radius'], 
        allow_overlap=params['allow_overlap'],
        verbose=params['verbose'], 
        standardize=params['standardize'], 
        t_r=params['t_r'], 
        high_pass=params['high_pass'], 
        low_pass=params['low_pass'], 
        detrend=params['detrend'],
        smoothing_fwhm=params['smoothing_fwhm']
    )
    if regressor_file is not None:
        extractor.set_regressors(regressor_file, params['regressors'], 
                                 params["load_confounds_kwargs"])

    if (params['discard_scans'] is not None) and (params['discard_scans'] > 0):
        extractor.discard_scans(params['discard_scans'])
    
    extractor.extract()
    out = os.path.join(params['out_dir'], replace_file_ext(input_file))
    extractor.save(out, params['n_decimals'])

    return out, extractor

def main():
    params = vars(_cli_parser())
    params = _check_nifti_params(params)
    metadata_path = make_param_file(params)
    shutil.copy2(params['roi_file'], metadata_path)

    # setup and run extraction
    run_extraction(extract_nifti, params['input_files'], params['roi_file'], 
                   params)

if __name__ == '__main__':
    raise RuntimeError("`nixtract/cli/nifti.py` should not be run directly. "
                       "Please `pip install` nixtract and use the "
                       "`nixtract-nifti` command.")
