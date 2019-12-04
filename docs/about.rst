
.. _outputs: outputs.html
__ outputs_
.. image:: ../resources/logo.png

`Niimasker` is a simple command-line wrapper for `Nilearn`'s `Masker objects <https://nilearn.github.io/manipulating_images/masker_objects.html>`_ (hence the name), which let you easily extract out region-of-interest (ROI) timeseries from functional MRI data while providing several options for applying additional post-processing (e.g., spatial smoothing, temporal filtering, confound regression, etc). This tool ultimately aims to extend many of `nilearn`'s powerful and convenient masking features to non-Python users (e.g., R and MATLAB users) who wish to analyze fMRI data.

In addition to providing a CLI for `Nilearn`'s masking functions, `Niimasker` also generates reports that document the extraction process. A report is generated for each functional image that `Niimasker` processes, and several plots are displayed for quality inspection. The goal is to provide a completely transparent account of your data extraction and post-processing step in your fMRI pipeline, while providing a simple and intuitive interface. To learn more about the outputs of Niimasker, see `Outputs of Niimasker`__.  

If you are using this project, please cite `Nilearn`:

Abraham, A., Pedregosa, F., Eickenberg, M., Gervais, P., Mueller, A., Kossaifi, J., … Varoquaux, G. (2014). `Machine learning for neuroimaging with scikit-learn. <https://www.frontiersin.org/articles/10.3389/fninf.2014.00014/full>`_ *Frontiers in Neuroinformatics*, *8*, 14.

See also: `Nilearn documentation <https://nilearn.github.io/index.html>`_.