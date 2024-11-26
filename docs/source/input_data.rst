########################
User Input File
########################

.. note:: A future iteration of the user input file will change to an excel file to make it easier. It is currently a text file.


The stellar input file contains all the metadata required for each star whose spectral energy distribution (SED) you wish to model. In the JSON parameterfile (file github | page description xxx links), the file path is described as

.. code-block:: rst

    ``['folders']['topdir']``/``['folders']['supportdir']``/StellarInputFiles

where, in the parameter file, and in the online github example, it's

::

    twinkle-master/Inputs_and_Models/StellarInputFiles/sample_stardata.txt.

The input file is a tab separated value (tsv) file, with a header where each column corresponds to a different meta-data parameter for each star. Each subsequent line contains information for a different star.

Below is a simple stellar input file (broken up for easier viewing, but it would be a 3 row by 20 column file) that includes optical photometric data and some mid-infrared `WISE data <https://www.jpl.nasa.gov/missions/wide-field-infrared-survey-explorer-wise/>`_. You can download The contents of this simple input file here (xxx github link).


   +-----------+-----+--------+----------+--------+------+
   | MainName  | spt | SPT2HIP| NoOptical| model  | temp |
   +===========+=====+========+==========+========+======+
   | Star_ID1  | A2V | A      | FALSE    | NextGen| 8000 |
   +-----------+-----+--------+----------+--------+------+
   | Star_ID2  | B5V | F      | TRUE     | ATLAS9 | 8000 |
   +-----------+-----+--------+----------+--------+------+

   +-----+-----+-----+-------+-----+------+-----+------+
   | grav| met | plx | e_Plx | B-V | e_B-V| BJm | BJme |
   +=====+=====+=====+=======+=====+======+=====+======+
   | 40  | 0   | 15.6| 0.25  | 0.07| 0.001| 5.2 | 0.01 |
   +-----+-----+-----+-------+-----+------+-----+------+
   | 45  | 0   | 20  | 0.3   | 0.1 | 0.001| 5.7 | 0.01 |
   +-----+-----+-----+-------+-----+------+-----+------+

   +-----+------+----+------+-----+------+
   | VJm | VJme | W1m| W1me | W2m | W2me |
   +=====+======+====+======+=====+======+
   | 5.1 | 0.01 | 4.9| 0.07 | 4.6 | 0.04 |
   +-----+------+----+------+-----+------+
   | 5.9 | 0.01 | 6.2| 0.07 | 4.6 | 0.04 |
   +-----+------+----+------+-----+------+

The file, shows meta-data for two stars, `Star_ID1`, and `Star_ID2` under the ``MainName`` column. All the columns up to ``e_B-V`` are not optional. This of course implies that the photometric columns are optional. You don't technically need to include any photometric bands, besides Johnson `B` & `V` but then you wouldn't get anything modeled. So you SHOULD include SOME photometric bands.

.. important::
    For the moment, one must include the Johnson B and V photometry and uncertainties (``BJm``, ``BJme``, ``VJm``, ``VJme``).


A full stellar input file can be found at this github link (xxx insert link). You'll see a number of different columns that you can use for your simulation. The file ``Input_StarFile_Description.xlsx`` contains information for all acceptable meta-data columns in your input stellar file, as well as column `descriptions`, `data type`, `units`, `parameter restrictions`, `examples` of the data, and whether the column is optional for the simulation or not. A copy of the file can be found below and at this link (xxx insert github link).

.. raw:: html

   <div style="text-align: center;">
       Description of meta-data that can be found and used in the stellar input file.<br>
       Scroll vertically and horizontally to see the full file.
   </div>

.. raw:: html

   <div style="max-height: 400px; overflow-y: scroll; overflow-x: scroll; border: 2px solid #ccc; padding: 15px;">


.. xlsx-table::
    :file: ../../Inputs_and_Models/StellarInputFiles/Input_StarFile_Description.xlsx
    :start-row: 2
    :start-column: 1
    :header-rows: 1
    :include-rows: 2-59
    :include-columns: 1-8


.. raw:: html

   </div>


The photometric bands that can be included for SED and excess flux calculations are listed in the  file under ``/Inputs_and_Models/RSR/available_filters.txt`` (XXX insert github link to file).

To see what the output of the modeling would look like with different meta-data in the stellar input file, check out the Jupyter Notebook tutorial at (XXX insert link).

.. important::
    The column strings should not have the asterisks. If the * columns are not included,
    then the ** columns in the "optional" column are required. If the *** columns are included, then "changekeys" in the JSON file must be set to "true" for these columns to be used.

.. important::
    To include spectral data columns, the parameter names should be in the following format: ``[band]m``, ``[band]me``, ``[band]_flux``, ``[band]_fluxe``. If the data is photometric, use ``[band]m``, and ``[band]me``, and the other two for fluxes in Jy.


