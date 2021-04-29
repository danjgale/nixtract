
import os
import argparse
import natsort

from nixtract.qc.analysis import quality_analysis
from nixtract.qc.utils import read_dist_file
# from .report import make_report


def _cli_parser():
    """Reads command line arguments and returns input specifications"""
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group(title='required arguments')
    required.add_argument('-t', "--timeseries", nargs="+", type=str, required=True,
                          help="One or more tab-delimited timeseries files. Can "
                               "also be a single string with a wildcard (*) to "
                               "specify all files matching the file pattern. If "
                               "so, these files are naturally sorted by file name "
                               "prior to extraction. Timeseries files must have "
                               "column headers in the first row")
    required.add_argument("-c", "--confounds", nargs="+", type=str, required=True,
                          help="One or more tab-delimited confounds file. Can "
                               "also be a single string with a wildcard (*) to "
                               "specify all files matching the file pattern. If "
                               "so, these files are naturally sorted by file name "
                               "prior to extraction. Each column should contain a "
                               "confound timeseries, and the first row should be "
                               "the confound name (i.e. column header). REQUIRED "
                               "CONFOUNDS: 6 head motion parameters: 'trans_x', "
                               "'trans_y' 'trans_z', 'rot_x', 'rot_y', 'rot_z'; "
                               "Framewise displacement 'framewise_displacement'; "
                               "DVARS `dvars`")
    required.add_argument("-o", "--out_dir", type=str, required=True,
                          help="The path to the output directory. Created if it"
                               "does not already exist")
    parser.add_argument("-d", "--distance_file", type=str,
                        help="File used to compute distances between regions "
                             "in the timeseries. Can either be 1) a volumetric "
                             "NIfTI image (.nii.gz or .nii) or 2) a "
                             "tab-delimited (.tsv) coordinate file. If a NIfTI "
                             "image, then region labels should be in the same "
                             "order as the timeseries (i.e. region 1 should be "
                             "column 1). If a .tsv file, must include columns "
                             "'x', 'y', 'z' (other columns will be ignored), "
                             "with each region as a row. Rows must be in the "
                             "same order as the timeseries file (i.e. row 1 "
                             "must be column 1 of the timeseries)")
    parser.add_argument('--group_only', type=int, default=1,
                        help='Skip plots for individual timeseries and only '
                             'generate group-level plots')   
    parser.add_argument('--n_jobs', type=int, default=1,
                        help='The number of CPUs to use if parallelization is '
                             'desired. Default: 1 (serial processing)')
    parser.add_argument('-v', '--verbose', action='store_true', default=False, 
                        help='Print out extraction progress')                  
    return parser.parse_args()


def _check_files(files):
    if len(files) == 1:
        # either a single file or cli glob did not expand
        if not os.path.exists(files[0]):
            raise ValueError('Provided timeseries file(s) do not exist. Check '
                             '-t (--timeseries)')
    return natsort.natsorted(files)


def main():

    params = vars(_cli_parser())
    timeseries = _check_files(params['timeseries'])
    confounds = _check_files(params['confounds'])

    if params['distance_file']:
        coordinates = read_dist_file(params['distance_file'])
    else:
        coordinates = None

    quality_analysis(timeseries, confounds, coordinates, params['out_dir'], 
                     n_jobs=params['n_jobs'], verbose=params['verbose'])
    # make_report(params['out_dir'])