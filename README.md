# deblend_sofia_detections

=====

Introduction
------------
    A  package to deblend Sofia detections based on the method described in 
    https://ui.adsabs.harvard.edu/abs/2025ApJ...980..157H/abstract

    Most of this is an adapted copy of the accompanying ipython notebook of that paper 
    and as such Qifeng Huang should be considered as a main author of this code


Requirements
------------
The code requires full installation of:

    python v3.11 or higher
    SoFiA-2
    
[python](https://www.python.org/)
[SoFiA-2](https://gitlab.com/SoFiA-Admin/SoFiA-2)


Installation
------------

Download the source code from the Github or simply install with pip as:

  	pip install deblend_sofia_detections

This should also install all required python dependencies.
We recommend the use of python virtual environments. If so desired a deblend_sofia_detections installation would look like:

  	python3 -m venv deblend_sofia_detections_venv

  	source deblend_sofia_detections_venv/bin/activate.csh

    pip install deblend_sofia_detections
(In case of bash the correct middle line is source deblend_sofia_detections_venv/bin/activate)

You might have to re-source the env:

  	source deblend_sofia_detections_venv/bin/activate.csh

Once you have installed deblend_sofia_detections you can check that it has been installed properly by running deblend_sofia_detections as.

  	deblend -v 


Running deblend_sofia_detections_venv
------------------

You can run deblend_sofia_detectionsby providing a configuration file by 

    deblend configuration_file=file.yml

an example yaml file with all parameters can be printed by running

    deblend print_examples=true 

If you only have a mask and a cube you can deblend those too with:

    deblend_cube sofia.original_data_cube=<HI data cube.fits> sofia.original_mask=<blended mask.fits>

please see the advanced input in readthedocs for an explanation of all parameters.