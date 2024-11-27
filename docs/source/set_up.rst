Set-Up
===========

.. toctree::
   :maxdepth: 2
   :caption: Contents:


JSON Parameterfile
*******************

The heart of `Twinkle` lies in the input parameter file, which is in JSON format. This file is essential for configuring `Twinkle`, including specifying the data folder structure, the names of relevant data files, the photometric bands used to fit the stellar spectral energy distribution (SED), plotting options, and more.

.. _json_parameterfile_link-label:

You can find an example JSON parameter file for `Twinkle` here: `JSON Parameter File <https://github.com/astropatel/twinkle/blob/master/paramfile.json>`_.

The JSON parameter file is organized as a dictionary. It contains boolean values, as well as dictionaries whose members are booleans, strings, arrays, integers, or floats. Below is the full JSON file, followed by descriptions of each section:


.. raw:: html

   <div style="max-height: 400px; overflow-y: scroll; border: 1px solid #ccc; padding: 10px;">

.. literalinclude:: paramfile.json
   :language: json
   :caption: JSON Parameterfile for Twinkle

.. raw:: html

   </div>

.. raw:: html

   <div style="text-align: center;">
       JSON Input Parameterfile for Twinkle
   </div>



comments
----------

The **comments** section contains metadata descriptions for various elements in the parameter file. These also describe the elements that are not under the top-level keys: ``longwave_bool``, ``W3Adapt``, ``W2Adapt``, ``fitphot``, ``scalephot``, ``changekeys``, ``satcheck``.



folders | files
--------------------

This section contains the names of relevant directories and files for running the code, as well as booleans for writing data to file. In parantheses are the current names listed in the example parameter file.

**Directory names**:

- ``topdir``: The top-level directory for the data (`twinkle-master`).
- ``supportdir``: Directory for support files (`Inputs_and_Models`).

**Empirical File Names**:

- ``starfile``: Filename for stellar data (`sample_stardata.txt`).
- ``input_stars``: Filename for input stars (`select_stars.txt`).
- ``bv_colorfile``: Empirical color file (`EMamajek_MSColors.txt`).

The names of all the parameters so far can be changed to whatever your heart desires.

**Write booleans parameters**:

- ``write2file``: Boolean indicating whether to write to a file.
- ``write_sed``: Boolean indicating whether to write SED data to a file.


plot
---------

The **plot** section contains options for visualizing the stellar data and SED.

- ``plot_any``: Boolean to enable or disable plotting.
- ``plot_grid``: Boolean to enable grid plots.
- ``plot_single``: Boolean to enable single plots.
- ``pRow`` and ``pCol``: Grid dimensions for the plot (rows and columns - ints).
- ``plotsize``: Dimensions of the plot (width, height - floats).
- ``ylim_up`` and ``ylim_low``: The upper and lower limits for the plot’s y-axis.
- ``xmax``: Maximum value for the x-axis.

spec_sample
--------------------

The **spec_sample** section defines parameters for sampling the spectrum.

- ``wave_min``: Minimum wavelength in angstroms (e.g., 2000Å - float).
- ``wave_max``: Maximum wavelength for calculations (Å - float).
- ``gridpts``: Number of grid points to sample across the spectrum (int).

phot
-----------

The **phot** section contains photometric band arrays used for various calculations. as well as two boolean parameters.

Wavelength array parameters:

- ``mags2use0_original``: Photometry bands used in the session.
- ``mags4Phot0_original``: Photometry bands used for photospheric fits.
- ``mags4scale0_original``: Bands used for scaling the raw photospheric spectrum for guessing the initial fit.
- ``mags4Dust0``: Bands used to calculate blackbody from excess flux.
- ``mag4DustWFC``: Bands for WFC dust calculations.
- ``scaleSEDbands``: Bands used to scale the final SED.
- ``Remove_RedStars``: Bands to not include for late-type stars.

Each of these is an array of strings of photometric bands that the code will use in different scenarios.

For example, if I set:

.. code-block:: json

    {
     "mags2use0_original" : ["BJ","VJ","J2M","H2M","Ks2M","W1","W2","W3","W4","HPACS160_flux","HPACS100_flux"],
    }

Then the code will allow the use of the Johnson **B**, and **V** bands, 2MASS **J**, **H**, **Ks**, and WISE All-Sky bands 1 through 4 in the simulation.

If you set ``mags4Phot0_original`` to only use the optical and short-wave IR bands, then those are the bands that will be used to fit the photospheric data.

.. important::
    You can insert any string combination of band descriptors into the different ``phot`` parameters in the parameterfile, so long as it also exists in ``mags2use0_original``.

.. important::
    The photometric descriptors used in the string arrays can be found in the `/Inputs_and_Models/RSR/Filters_README.txt <https://github.com/astropatel/twinkle/blob/master/Inputs_and_Models/RSR/Filters_ReadME.txt>`_  and is listed below. This file also contains information on the data sources for these relative spectral response curves. Captilization matters of each file matters.

.. _rsr_filters_readme-label:

.. raw:: html

   <div style="max-height: 400px; overflow-y: scroll; border: 1px solid #ccc; padding: 10px;">

.. literalinclude:: _static/Filters_ReadME.txt
   :language: text
   :caption: Supported Spectral Response Filters

.. raw:: html

   </div>

.. raw:: html

   <div style="text-align: center;">
       Available Photometric Filters User Can Use
   </div>



WISE_excess
-------------

The **WISE_excess** section defines WISE photometry cuts for identifying excess in infrared bands:

- ``W13_cut``: Cutoff for W1 and W3 excess.
- ``W23_cut``: Cutoff for W2 and W3 excess.
- ``W12_cut``: Cutoff for W1 an d W2 excess.

Data Files and Folder Structure
******************************************

.. _directory_structure-label:

Directory Structure
---------------------

Once you've installed `Twinkle` and it is included in your path, you'll need to set-up your working-folder which will contain input data to model the stellar spectral energy distribution of a given target, and the empirical data the code will use to perform that calculation.

You can download the current data and folders from the `Github project page <https://github.com/astropatel/twinkle/tree/master>`_. The `twinkle` package need not be in your working directory as long as it's in your ``PYTHONPATH``.

The data should be in the structure shown in the first panel below. The code will use information from the ``folders`` and ``files`` elements of the :ref:`JSON parameterfile <json_parameterfile_link-label>` to identify which folder to pull data from. The color-coding in the below directory structure is as follows:

.. raw:: html

    <ul>
        <li><span style="color: red;">Folder name from the parameter file</span></li>
        <li><span style="color: lightcoral;">Static folder names</span></li>
        <li><span style="color: blue;">File name from the parameter file</span></li>
    </ul>

You need to make sure that ``folders``/``topdir`` points to where the folder is on your machine.

 todo: make paramfile auto-read from directory structure.

.. raw:: html

   <div style="text-align: center;">
       Twinkle Directory Structure using Parameterfile placeholders
   </div>

.. raw:: html

   <div style="max-height: 400px; overflow-y: scroll; border: 1px solid #ccc; padding: 10px;">

    <pre style="font-size:70%;">
    <span style="color: red;">['folders']['topdir']</span>
        ├── Twinkle Tutorial.ipynb
        ├── paramfile.json
        ├── <span style="color: red;">['folders']['supportdir']</span>
        │   ├── <span style="color: blue;">['files']['bv_colorfile']</span>
        │   ├── <span style="color: lightcoral;">StellarInputFiles/</span>
        │   │   ├── Input_StarFile_Description.xlsx
        │   │   └── <span style="color: blue;">['files']['starfile']</span>
        │   ├── <span style="color: lightcoral;">StellarGridModels/</span>
        │   │   ├── <span style="color: lightcoral;">ATLAS9/</span>
        │   │   │   ├── README.txt
        │   │   │   ├── <span style="color: lightcoral;">km05/</span>
        │   │   │   │   ├── km05_10000.fits
        │   │   │   │   ├── km05_10500.fits
        │   │   │   │   .
        │   │   │   │   .
        │   │   │   │   .
        │   │   │   │   ├── km05_9500.fits
        │   │   │   │   └── km05_9750.fits
        │   │   │   └── <span style="color: lightcoral;">kp00/</span>
        │   │   │       ├── kp00_10000.fits
        │   │   │       ├── kp00_10500.fits
        │   │   │       .
        │   │   │       .
        │   │   │       .
        │   │   │       ├── kp00_9500.fits
        │   │   │       └── kp00_9750.fits
        │   │   └── <span style="color: lightcoral;">NextGen/</span>
        │   │       ├── lteNextGen_10000.0_3.5_0.0.txt
        │   │       ├── lteNextGen_10000.0_4.0_0.0.txt
        │   │       ├── lteNextGen_10000.0_4.5_0.0.txt
        │   │       .
        │   │       .
        │   │       .
        │   │       └── lteNextGen_9800.0_5.5_0.0.txt
        │   │  
        └───└── <span style="color: lightcoral;">RSR/</span>
                ├── Akari90.dat
                ├── Akari_IRCL18W.dat
                ├── Akari_IRCS9W.dat
                 .
                 .
                 .
                ├── W1_WISE.dat
                ├── W2_WISE.dat
                ├── W3_WISE.dat
                ├── W4_WISE.dat
                └── Filters_ReadME.txt
    </pre>

.. raw:: html

   </div>

.. raw:: html

   <div style="text-align: center;">
       <br>
   </div>


For example, if we take the values of the ``files`` and ``folders`` dictionaries in the example parameter file on this page, the directory structure will look like the image below. Now, you might ask, why is the **Inputs_and_Models** folder under the ``folders``/``supportdir`` a variable? Why would you need to change that name.

Sure. You could ask that.

** can topdir be made to be a a file path?**

.. raw:: html


   <div style="text-align: center;">
       Twinkle Directory Structure using Parameterfile values
   </div>

.. raw:: html

   <div style="max-height: 400px; overflow-y: scroll; border: 1px solid #ccc; padding: 10px;">

    <pre style="font-size:70%;">
    <span style="color: red;">twinkle-master</span>
        ├── Twinkle Tutorial.ipynb
        ├── paramfile.json
        ├── <span style="color: red;">Inputs_and_Models/</span>
        │   ├── <span style="color: blue;">EMamajek_MSColors.txt</span>
        │   ├── <span style="color: lightcoral;">StellarInputFiles/</span>
        │   │   ├── Input_StarFile_Description.xlsx
        │   │   └── <span style="color: blue;">sample_stardata.txt</span>
        │   ├── <span style="color: lightcoral;">StellarGridModels/</span>
        │   │   ├── <span style="color: lightcoral;">ATLAS9/</span>
        │   │   │   ├── README.txt
        │   │   │   ├── <span style="color: lightcoral;">km05/</span>
        │   │   │   │   ├── km05_10000.fits
        │   │   │   │   ├── km05_10500.fits
        │   │   │   │   .
        │   │   │   │   .
        │   │   │   │   .
        │   │   │   │   ├── km05_9500.fits
        │   │   │   │   └── km05_9750.fits
        │   │   │   └── <span style="color: lightcoral;">kp00/</span>
        │   │   │       ├── kp00_10000.fits
        │   │   │       ├── kp00_10500.fits
        │   │   │       .
        │   │   │       .
        │   │   │       .
        │   │   │       ├── kp00_9500.fits
        │   │   │       └── kp00_9750.fits
        │   │   └── <span style="color: lightcoral;">NextGen/</span>
        │   │       ├── lteNextGen_10000.0_3.5_0.0.txt
        │   │       ├── lteNextGen_10000.0_4.0_0.0.txt
        │   │       ├── lteNextGen_10000.0_4.5_0.0.txt
        │   │       .
        │   │       .
        │   │       .
        │   │       └── lteNextGen_9800.0_5.5_0.0.txt
        │   │  
        └───└── <span style="color: lightcoral;">RSR/</span>
                ├── Akari90.dat
                ├── Akari_IRCL18W.dat
                ├── Akari_IRCS9W.dat
                 .
                 .
                 .
                ├── W1_WISE.dat
                ├── W2_WISE.dat
                ├── W3_WISE.dat
                ├── W4_WISE.dat
                └── Filters_ReadME.txt
    </pre>


.. raw:: html

   </div>

.. raw:: html

   <div style="text-align: center;">
       <br>
   </div>


['folders']['topdir'] | (e.g., twinkle-master/)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Top-level directory that contains all the input and empirical data Twinkle will use for the simulation.

Twinkle Tutorial.ipynb
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Jupyter notebook tutorial that will help you get started.


['folders']['supportdir'] | (e.g., Inputs_and_Models/)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This folder contains the synthetic stellar grid models, relative spectral response (RSR) files, empirical stellar color relations, and descriptions of input data.

RSR/
^^^^^

Relative spectral response (RSR) data corresponding to the identifiers in the available photometric filters file (link to table above).

* :ref:`Available filters list<rsr_filters_readme-label>`
* `Data for all the filters <https://github.com/astropatel/twinkle/tree/master/Inputs_and_Models/RSR>`_
* :ref:`Description for RSR folder and filters<RSR_Description-label>`

['files']['bv_colorfile'] | (e.g., EMamajek_MSColors.txt)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This file contains empirical color relations for all stellar spectral types. The data is taken from `Dr. Eric Mamajek's <https://www.pas.rochester.edu/~emamajek/>`_ carefully maintained `color tables <https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt>`_.

The file used by `Twinkle` can be found `on the github page <https://github.com/astropatel/twinkle/blob/master/Inputs_and_Models/EMamajek_MSColors.txt>`_.

.. warning:: Dr. Eric Mamajek's tables are constantly being updated, and the table used by the simulation may be outdated. Take this into consideration when using the table on Github.

['files']['starfile'] | (e.g., sample_stardata.txt, sample_stardata_simple.txt)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This file contains tabular tsv data for each star whose SED you wish to model. Each row contains meta-data and empirical photometric information for a single star, and the columns correspond to different meta-data information. Information on how to build the starfile can be found on the :doc:`User Input File page<input_data>`.

.. note:: future version of the user input file will be an Excel workbook.

Input_StarFile_Description.xlsx
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

File containing descriptions of the meta data for the stellar input data the code will use in the grid model fitting. `Here's the link <https://github.com/astropatel/twinkle/blob/master/Inputs_and_Models/~%24Input_StarFile_Description.xlsx>`_.

You can also check out the contents on the User Input File page :ref:`here<input_file_description_image-label>`

StellarGridModels/
^^^^^^^^^^^^^^^^^^^^^^^^

Contains atmospheric stellar grid models. Links to the models (ATLAS9 and NextGen) on the Github page are `here <https://github.com/astropatel/twinkle/tree/master/Inputs_and_Models/StellarGridModels>`_.

Description to the atmopsheric models can be found on the :doc:`Model Data<model_data>` page.

ATLAS9/
^^^^^^^^

These are the `Kurucz Atlas 9 <https://ui.adsabs.harvard.edu/abs/1993yCat.6039....0K/abstract>`_ models and are in FITS format. More information can be found in the :ref:`ATLAS9 section on the Model Data page<atlas9-label>`.

NextGen/
^^^^^^^^^

These are the `NextGen` atmospheric models (`Hauschildt et al. 1999 <https://iopscience.iop.org/article/10.1086/306745>`_) models. More information can be found in the :ref:`NexstGen section on the Model Data page<nextgen-label>`.






