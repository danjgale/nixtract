=====
Usage
=====

`Niimasker` is run using the command-line, but is flexibly parameterized using either the command-line itself, or a JSON configuration file. All parameters except ``output_dir`` can be either specified as command-line arguments or as parameters in the configuration file. The parameters set in the configuration file take priority over those called via the command-line(i.e. parameters set in the configuration file will overwrite those set in the CLI). **The configuration file is the preferred way of using `niimasker`**; refer to `The Configuration File`_ for more detail.  

Command-Line Arguments
======================

.. argparse::
   :ref: niimasker.cli._cli_parser
   :prog: niimasker

Many of the parameters map directly onto the masker parameters in `nilearn` (see the `documentation <https://nilearn.github.io/modules/reference.html#module-nilearn.input_data>`_ and `user guide <https://nilearn.github.io/building_blocks/manual_pipeline.html#masking>`_ for more detail). The underlying masker object is determined by the `--as_voxels` flag, which will set the masker object to ``nilearn.input_data.NiftiMasker`` if ``True``. Otherwise, ``nilearn.input_data.NiftiMapsMasker`` is used to obtain averaged region timecourses.  

**Required arguments**
-  ``ouput_dir``, specified by command-line only
- ``input_files``, can be specified by the command-line or by a configuration file
- ``mask_img``, can be specified by the command-line or by a configuration file

All other arguments are optional.


The Configuration File
=======================

Instead of passing all of the parameters through the command-line, `niimasker` also provides support for a simple configuration JSON file. An empty configuration file template of all of the parameters is shown below:

.. code-block:: none

  {
    "input_files": [],
    "mask_img": "",
    "labels": [],
    "regressor_files": null,
    "regressor_names": [],
    "as_voxels": false,
    "standardize": false,
    "t_r": null,
    "detrend": false,
    "high_pass": null,
    "low_pass": null,
    "smoothing_fwhm": null,
    "discard_scans": null,
    "n_jobs": 1
  }


