# Package Name

=====

Introduction
------------



Requirements
------------
The code requires full installation of:

    python v3.6 or higher
    


[python](https://www.python.org/)


Installation
------------

Download the source code from the Github or simply install with pip as:

  	pip install packagename

This should also install all required python dependencies.
We recommend the use of python virtual environments. If so desired a TRM_errors installation would look like:

  	python3 -m venv package_name_venv

  	source package_name_venv/bin/activate.csh

    pip install package_name

(In case of bash the correct middle line is 	source package_name_venv/bin/activate)
You might have to restart the env:

  	deactivate

  	source package_name_venv/bin/activate.csh

Once you have installed FAT you can check that it has been installed properly by running FAT as.

  	executable -v 


Running package_name_venv
------------------

You can run package_name by providing a configuration file by 

executable configuration_file=file.yml

an example yaml file with all parameters can be printed by running

executable print_examples=true 

please see the advanced input in readthe docs for an explanation of all parameters.