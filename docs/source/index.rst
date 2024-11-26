.. twinkle documentation master file, created by
   sphinx-quickstart on Tue Oct 15 22:19:40 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

=========
TWINKLE! 
=========


Introduction
*************
*Twinkle* is a Python-based tool that has the ability to calculate the spectral energy distribution (`SED <https://coolwiki.ipac.caltech.edu/index.php/SED_plots_introduction>`_) of a stellar source using empirical photometric data and stellar model grids.

*Twinkle* was originally created to help me calculate the excess infrared (IR) flux from a star. The presence of an IR excess is evidence for dust orbiting the star. This dust would have been created by the grinding and collisions of asteroids, usually moved around through the gravitational influence of a larger planetary object -- so basically the presence of an IR excess points to the potential for finding planets. You can check out the published papers from my thesis using this code in `Patel, Metchev, and, Heinze, 2014 <https://iopscience.iop.org/article/10.1088/0067-0049/212/1/10>`_ and `Patel, et al., 2017 <https://iopscience.iop.org/article/10.3847/1538-3881/153/2/54>`_.

Interested in learning more about debris disks? Check out my `blog post <http://cosmicdiary.org/geminiplanetimager/2015/03/04/debris-disks-searching-for-dust-to-find-planets/>`_.

This code base can help you quickly calculate the temperature and location of the dust to first order by fitting the assumed blackbody or modified blackbody function to the broadband excess emission.

Feel free to fork and contribute.


Available Features
*******************

- Model multiple stellar SEDs at once with easy to use stellar input file.

- Late B to K-spectral type modeling support. Additional spectal type modeling is possible but not tested (*yet*)

- Plotting capabilities (link to plotting page upcoming) of the empirical data, and the modeled distribution.

- This `Jupyter Notebook <https://github.com/astropatel/twinkle/blob/master/Twinkle_Tutorial.ipynb>` will give you a quick-start on using `Twinkle`.


.. note:: The code hasn't been tested in a while so there will probably be issues with it that are still being worked on. If you want to contribute, fork and pull-request things up.


.. toctree::
   :numbered:
   :maxdepth: 1
   :caption: Contents:

   set_up
   input_data
   model_data
   modules



.. :class:`sed.blahblah`
.. :meth:`sed.blahblah`


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
