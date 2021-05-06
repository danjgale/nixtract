# nixtract

<p align="center">
  <img src="resources/logo.png" alt="logo" width="500"/>
</p>

`nixtract` (**N**euro**I**maging e**XTRACT**ion) is a collection of simple command-line tools that provide a decently unified interface to extract and process timeseries data from NIFTI, GIFTI, and CIFTI neuroimaging files. 

The primary goal of `nixtract` is to provide the user with ready-to-use timeseries data for a variety of fMRI analyses. `nixtract` can extract the mean timeseries for each region in a provided atlas, or `nixtract` can also extract the timeseries of individual voxels/vertices within a specified region. These timeseries can be flexibly denoised using temporal filtering/detrending, spatial smoothing, and confound regression, thus providing the user fully processed timeseries for subsequent analysis.   

## Installation

`nixtract` can be installed via `pip` using the following command:

```
pip install nixtract
```

## Command Line Interfaces

Nixtract has a command-line interface (CLI) for each file type, as overviewed below:

### NIFTIs

Data can be extract from NIFTI (`.nii` or `.nii.gz`) data using `nixtract-nifti`:

```
usage: nixtract-nifti [-h] [--input_files INPUT_FILES [INPUT_FILES ...]]
                      [--roi_file ROI_FILE] [--mask_img MASK_IMG]
                      [--labels LABELS [LABELS ...]] [--as_voxels]
                      [--radius radius] [--allow_overlap]
                      [--smoothing_fwhm SMOOTHING_FWHM]
                      [--regressor_files REGRESSOR_FILES [REGRESSOR_FILES ...]]
                      [--regressors REGRESSORS [REGRESSORS ...]]
                      [--standardize] [--t_r T_R] [--high_pass HIGH_PASS]
                      [--low_pass LOW_PASS] [--detrend]
                      [--discard_scans DISCARD_SCANS] [--n_jobs N_JOBS]
                      [--n_decimals N_DECIMALS] [-c CONFIG] [-v]
                      out_dir

positional arguments:
  out_dir               The path to the output directory. Created if itdoes
                        not already exist

optional arguments:
  -h, --help            show this help message and exit
  --input_files INPUT_FILES [INPUT_FILES ...]
                        One or more input NIFTI images (.nii.gz only). Can
                        also be a single string with wildcards (*) to specify
                        all files matching the file pattern. If so, these
                        files are naturally sorted by file name prior to
                        extraction
  --roi_file ROI_FILE   Parameter that defines the region(s) of interest. This
                        can be 1) a file path to NIFTI image that is an atlas
                        of multiple regions, a probabilistic atlas, or a
                        binary mask of one region, 2) a nilearn query string
                        formatted as `nilearn:<atlas-name>:<atlas-parameters>
                        3) a file path to a .tsv file that has x, y, z columns
                        that contain coordinates in MNI space. Refer to online
                        documentation for more detail and how these options
                        map onto the underlying nilearn masker classes
  --mask_img MASK_IMG   File path of a binary mask a to be used when
                        `roi_file` is a) an multi-region atlas or b) a list of
                        coordinates. This will restrict extraction to only
                        voxels within the mask. If `roi_file` is a single
                        region binary mask, this will be ignored
  --labels LABELS [LABELS ...]
                        Labels corresponding to the region numbers in
                        `roi_file`. Can either be a) a list of strings, b) or
                        a .tsv file that contains a `Labels` column. Labels
                        must be sorted in ascending order to correctly
                        correspond to the atlas indices. The number of labels
                        provided must match the number of non-zero indices in
                        `roi_file`. Numeric indices are used if not provided
  --as_voxels           Extract the timeseries of each voxel instead in a a
                        region rather than the mean timeseries. This is only
                        available for single region (binary) masks. Default:
                        False
  --radius radius       Set the radius of the spheres (in mm) centered on the
                        coordinates provided in `roi_file`. Only applicable
                        when a coordinate .tsv file is passed to `roi_file`;
                        otherwise, this will be ignored. If not set, the
                        timeseries of each coordinate is extracted (nilearn
                        default)
  --allow_overlap       Permit overlapping spheres when coordinates are
                        provided to `roi_file` and `radius` is provided
  --smoothing_fwhm SMOOTHING_FWHM
                        Smoothing kernel FWHM (in mm) if spatial smoothing is
                        desired.
  --regressor_files REGRESSOR_FILES [REGRESSOR_FILES ...]
                        One or more tab-separated files with regressors in
                        each column. The number of files must match the number
                        of input files and must be in the same order. The
                        number of rows in each file must match the number of
                        timepoints in their respective input files. Can also
                        be a single string with a wildcard (*) to specify all
                        files matching the file pattern. If so, these files
                        are naturally sorted by file name prior to extraction.
                        Double check to make sure these are correctly aligned
                        with the input files (see the parameters.json in the
                        output)
  --regressors REGRESSORS [REGRESSORS ...]
                        Regressor names or strategy to use for confound
                        regression. Must be a) list of specified column names
                        in all of the regressor_files, b) a predefined
                        strategy by load_confounds, or c) a list compatible
                        with load_confounds flexible denoising strategy
                        options. See the documentation
                        https://github.com/SIMEXP/load_confounds. If no
                        regressor information provided but regressor files are
                        provided, then all regressors in regressor files are
                        used
  --standardize         Whether to standardize (z-score) each timeseries.
                        Default: False
  --t_r T_R             The TR of the functional files, specified in seconds.
                        Required if temporal filtering/detrending is specified
  --high_pass HIGH_PASS
                        High pass filter cut off in Hertz. Do not use if high
                        pass cosine regressors are specified in `regressors`
  --low_pass LOW_PASS   Low pass filter cut off in Hertz. Do not use if low
                        pass cosine regressors are specified in `regressors`
  --detrend             Temporally detrend the data. Default: False
  --discard_scans DISCARD_SCANS
                        Discard the first N scans of each functional image
  --n_jobs N_JOBS       The number of CPUs to use if parallelization is
                        desired. Default: 1 (serial processing)
  --n_decimals N_DECIMALS
                        Specify the number of decimals for output timeseries
                        files. Fewer decimals are recommended for reducing
                        disk-space, particularly for large extractions
  -c CONFIG, --config CONFIG
                        A configuration .json file to pass parameters This
                        will overwrite command-line arguments if the same
                        parameter is specified in both. See online
                        documentation for formatting and what keys to include
  -v, --verbose         Print out extraction progress
```

### GIFTIs

Data can be extract from GIFTI functional files (`.func.gii`) data using `nixtract-gifti`:

```
usage: nixtract-gifti [-h] [--lh_files LH_FILES [LH_FILES ...]]
                      [--rh_files RH_FILES [RH_FILES ...]]
                      [--lh_roi_file LH_ROI_FILE] [--rh_roi_file roi_file]
                      [--as_vertices] [--denoise-pre-extract]
                      [--regressor_files REGRESSOR_FILES [REGRESSOR_FILES ...]]
                      [--regressors REGRESSORS [REGRESSORS ...]]
                      [--standardize] [--t_r T_R] [--high_pass HIGH_PASS]
                      [--low_pass LOW_PASS] [--detrend]
                      [--discard_scans DISCARD_SCANS] [--n_jobs N_JOBS]
                      [--n_decimals N_DECIMALS] [-c CONFIG] [-v]
                      out_dir

positional arguments:
  out_dir               The path to the output directory. Created if itdoes
                        not already exist

optional arguments:
  -h, --help            show this help message and exit
  --lh_files LH_FILES [LH_FILES ...]
                        One or more input functional GIFTI files (.func.gii)
                        for the left hemisphere. Can also be a single string
                        with wildcards (*) to specify all files matching the
                        file pattern. If so, these files are naturally sorted
                        by file name prior to extraction
  --rh_files RH_FILES [RH_FILES ...]
                        One or more input functional GIFTI files (.func.gii)
                        for the right hemisphere. Can also be a single string
                        with wildcards (*) to specify all files matching the
                        file pattern. If so, these files are naturally sorted
                        by file name prior to extraction
  --lh_roi_file LH_ROI_FILE
                        A label GIFTI file (.label.gii) or a Freesurfer
                        annotation file (.annot) for the left hemipshere. Must
                        include one or more labels
  --rh_roi_file roi_file
                        A label GIFTI file (.label.gii) or a Freesurfer
                        annotation file (.annot) for the right hemipshere.
                        Must include one or more labels
  --as_vertices         Extract the timeseries of each vertex in a region
                        rather than the mean timeseries. This is only
                        available for when `lh_roi_file` and/or `rh_roi_file`
                        are a single region, i.e. a binary mask. Default:
                        False
  --denoise-pre-extract
                        Denoise data (e.g., filtering, confound regression)
                        before timeseries extraction. Otherwise, denoising is
                        done on the extracted timeseries, which is consistent
                        with nilearn and is more computationally efficient.
                        Default: False
  --regressor_files REGRESSOR_FILES [REGRESSOR_FILES ...]
                        One or more tab-separated files with regressors in
                        each column. The number of files must match the number
                        of input files and must be in the same order. The
                        number of rows in each file must match the number of
                        timepoints in their respective input files. Can also
                        be a single string with a wildcard (*) to specify all
                        files matching the file pattern. If so, these files
                        are naturally sorted by file name prior to extraction.
                        Double check to make sure these are correctly aligned
                        with the input files (see the parameters.json in the
                        output)
  --regressors REGRESSORS [REGRESSORS ...]
                        Regressor names or strategy to use for confound
                        regression. Must be a) list of specified column names
                        in all of the regressor_files, b) a predefined
                        strategy by load_confounds, or c) a list compatible
                        with load_confounds flexible denoising strategy
                        options. See the documentation
                        https://github.com/SIMEXP/load_confounds. If no
                        regressor information provided but regressor files are
                        provided, then all regressors in regressor files are
                        used
  --standardize         Whether to standardize (z-score) each timeseries.
                        Default: False
  --t_r T_R             The TR of the functional files, specified in seconds.
                        Required if temporal filtering/detrending is specified
  --high_pass HIGH_PASS
                        High pass filter cut off in Hertz. Do not use if high
                        pass cosine regressors are specified in `regressors`
  --low_pass LOW_PASS   Low pass filter cut off in Hertz. Do not use if low
                        pass cosine regressors are specified in `regressors`
  --detrend             Temporally detrend the data. Default: False
  --discard_scans DISCARD_SCANS
                        Discard the first N scans of each functional image
  --n_jobs N_JOBS       The number of CPUs to use if parallelization is
                        desired. Default: 1 (serial processing)
  --n_decimals N_DECIMALS
                        Specify the number of decimals for output timeseries
                        files. Fewer decimals are recommended for reducing
                        disk-space, particularly for large extractions
  -c CONFIG, --config CONFIG
                        A configuration .json file to pass parameters This
                        will overwrite command-line arguments if the same
                        parameter is specified in both. See online
                        documentation for formatting and what keys to include
  -v, --verbose         Print out extraction progress

```

### CIFTIs

Data can be extract from CIFTI functional files (`.dtseries.nii`) data using `nixtract-cifti`:

```
usage: nixtract-cifti [-h] [--input_files INPUT_FILES [INPUT_FILES ...]]
                      [--roi_file ROI_FILE] [--as_vertices]
                      [--denoise-pre-extract]
                      [--regressor_files REGRESSOR_FILES [REGRESSOR_FILES ...]]
                      [--regressors REGRESSORS [REGRESSORS ...]]
                      [--standardize] [--t_r T_R] [--high_pass HIGH_PASS]
                      [--low_pass LOW_PASS] [--detrend]
                      [--discard_scans DISCARD_SCANS] [--n_jobs N_JOBS]
                      [--n_decimals N_DECIMALS] [-c CONFIG] [-v]
                      out_dir

positional arguments:
  out_dir               The path to the output directory. Created if itdoes
                        not already exist

optional arguments:
  -h, --help            show this help message and exit
  --input_files INPUT_FILES [INPUT_FILES ...]
                        One or more input CIFTI dtseries files
                        (.dtseries.nii). Can also be a single string with
                        wildcards (*) to specify all files matching the file
                        pattern. If so, these files are naturally sorted by
                        file name prior to extraction
  --roi_file ROI_FILE   CIFTI dlabel file (.dlabel.nii) with one or more
                        labels
  --as_vertices         Extract the timeseries of each vertex in a a region
                        rather than the mean timeseries.This is only available
                        for when `roi_file` is single region, i.e. a binary
                        mask. Default: False
  --denoise-pre-extract
                        Denoise data (e.g., filtering, confound regression)
                        before timeseries extraction. Otherwise, denoising is
                        done on the extracted timeseries, which is consistent
                        with nilearn and is more computationally efficient.
                        Default: False
  --regressor_files REGRESSOR_FILES [REGRESSOR_FILES ...]
                        One or more tab-separated files with regressors in
                        each column. The number of files must match the number
                        of input files and must be in the same order. The
                        number of rows in each file must match the number of
                        timepoints in their respective input files. Can also
                        be a single string with a wildcard (*) to specify all
                        files matching the file pattern. If so, these files
                        are naturally sorted by file name prior to extraction.
                        Double check to make sure these are correctly aligned
                        with the input files (see the parameters.json in the
                        output)
  --regressors REGRESSORS [REGRESSORS ...]
                        Regressor names or strategy to use for confound
                        regression. Must be a) list of specified column names
                        in all of the regressor_files, b) a predefined
                        strategy by load_confounds, or c) a list compatible
                        with load_confounds flexible denoising strategy
                        options. See the documentation
                        https://github.com/SIMEXP/load_confounds. If no
                        regressor information provided but regressor files are
                        provided, then all regressors in regressor files are
                        used
  --standardize         Whether to standardize (z-score) each timeseries.
                        Default: False
  --t_r T_R             The TR of the functional files, specified in seconds.
                        Required if temporal filtering/detrending is specified
  --high_pass HIGH_PASS
                        High pass filter cut off in Hertz. Do not use if high
                        pass cosine regressors are specified in `regressors`
  --low_pass LOW_PASS   Low pass filter cut off in Hertz. Do not use if low
                        pass cosine regressors are specified in `regressors`
  --detrend             Temporally detrend the data. Default: False
  --discard_scans DISCARD_SCANS
                        Discard the first N scans of each functional image
  --n_jobs N_JOBS       The number of CPUs to use if parallelization is
                        desired. Default: 1 (serial processing)
  --n_decimals N_DECIMALS
                        Specify the number of decimals for output timeseries
                        files. Fewer decimals are recommended for reducing
                        disk-space, particularly for large extractions
  -c CONFIG, --config CONFIG
                        A configuration .json file to pass parameters This
                        will overwrite command-line arguments if the same
                        parameter is specified in both. See online
                        documentation for formatting and what keys to include
  -v, --verbose         Print out extraction progress
```

## The configuration JSON file

Instead of passing all of the parameters through the command-line, `nixtract` also provides support for a simple configuration JSON file. The only parameter that needs to be passed into the command-line is the output directory (`output_dir`). All other parameters can either be set by the configuration file or by the command-line. **Note that the configuration file overwrites any of the command-line parameters**. 

Not all parameters need to be included in the configuration file; only the ones you wish to use. An example use-case that combines both the command-line parameters and configuration file:

`nixtract-nifti output/ -i img_1.nii.gz img_2.nii.gz --config config.json`

Where `config.json` is:

```JSON
{
  "roi_file": "some_atlas.nii.gz",
  "standardize": true,
  "regressor_files": [
    "confounds1.tsv",
    "confounds2.tsv"
  ],
  "regressors": "Params6",
  "t_r": 2,
  "high_pass": 0.01,
  "smoothing_fwhm": 6
}
```
This set up is convenient when your `output_dir` and `input_files` vary on a subject-by-subject basis, but your post-processing and atlas might stay constant. Therefore, constants across subjects can be stored in the project's configuration file.

Configuration templates for each CLI are in `resources/config-templates`.  
