
=========
TWINKLE! 
=========


Introduction
*************
*Twinkle* is a Python-based tool designed to calculate the spectral energy distribution (`SED <https://coolwiki.ipac.caltech.edu/index.php/SED_plots_introduction>`_) of stellar sources using photometric data and stellar model grids.

Originally, I developed *Twinkle* to help me characterize excess infrared (IR) flux from stars. The presence of IR excess often indicates dust orbiting the star, likely produced by the grinding and collisions of asteroids, which are typically influenced by the gravity of larger planetary objects. In short, detecting IR excess can suggest the presence of planets. For more details, refer to the published papers from my thesis : `Patel, Metchev, and, Heinze, 2014 <https://iopscience.iop.org/article/10.1088/0067-0049/212/1/10>`_ and `Patel, et al., 2017 <https://iopscience.iop.org/article/10.3847/1538-3881/153/2/54>`_.

Interested in learning more about debris disks? Check out my `blog post <http://cosmicdiary.org/geminiplanetimager/2015/03/04/debris-disks-searching-for-dust-to-find-planets/>`_.

This codebase allows you to quickly estimate the temperature and location of the dust by fitting an assumed blackbody or modified blackbody function to broadband excess emission.

Feel free to fork and contribute! The core code and fitting routines were created over 10 years ago, so some functions might be outdated and could benefit from using modern ``SciPy`` alternatives. Updates are planned, but contributions are welcome—feel free to fix and push changes.





Available Features
*******************

- Model multiple stellar SEDs at once with easy to use stellar input file.

- Late B to K-spectral type modeling support. Additional spectal type modeling is possible but not tested (*yet*)

- Plotting capabilities (link to plotting page upcoming) of the empirical data, and the modeled distribution.

- This `Jupyter Notebook <https://github.com/astropatel/twinkle/blob/master/Twinkle_Tutorial.ipynb>`_ will give you a quick-start on using `Twinkle`.


.. note:: The code hasn't been tested in a while so there will probably be issues with it that are still being worked on. If you want to contribute, fork and pull-request things up.

Future Capabilities
---------------------

- Graphical User Interface

- Generation of several SEDs in one execution.

- Generation of (modified) blackbody fits to IR excess.


Installation & Quick-Start
****************************

Here is the `Twinkle` `Github Link <https://github.com/astropatel/twinkle>`_.

Feel free to fork it, install it, etc.

1. Install Twinkle
--------------------
You can install it through pip:

.. code-block:: bash

   pip install git+https://github.com/astropatel/twinkle.git

or by downloading it directly from the github project and adding it to your ``PYTHONPATH``.

2. Set-Up Directories and ParameterFile
-----------------------------------------

Next you'll need to make sure to get the directory structure of your working directory set-up right. This also includes editting the input parameterfile.

The parameterfile contains information on the names you've set for the  folders/files, photometric bands you will use to fit the SED, any associated excess, and more.

You can find directions on the set-up and parameterfile on the :doc:`Set-Up Page<set_up>`.

3. Create/Edit User Input File
-------------------------------

If you want your stars to `Twinkle`, you'll need to give the code information on which stars you want it process! The user input file is an excel file you can edit to include physical and photometric information on multiple stars the code can pull from.

For more information on the user input file, check out the :doc:`User Input File page<input_data>`.

4. Fire up the code
--------------------

Check out the tutorial from the `Jupyter Notebook <https://github.com/astropatel/twinkle/blob/master/Twinkle_Tutorial.ipynb>`_ on how to use `Twinkle` once you have everything set up.


.. toctree::
   :numbered:
   :maxdepth: 2
   :caption: Contents:

   set_up
   input_data
   model_data
   modules



The Logo
*********

I created the Twinkle logo after obtaining baseline ideas from AI image generators, and then drawing it on Inkscape. The top prong, or point, or (what is it called?) of the star image is meant to look like a stellar SED.


The Name
*********

I prefer my projects to have apt names, but then like to backronym it into something non-sensical and fun. While I haven't landed on one for `TWINKLE`, here are suggestions from friends and colleagues


* The WIse Non-gaseous disK Locator Extraordinaire (by: `Scott Barenfeld`)

* The Wise Infrared Non-gaseous disK Locator Extraordinaire (by: `Calen Henderson`)

* That's What I Need to Know about teLescopes from Experts (by: `Tiffany Meshkat`)



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
