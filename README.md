# Twinkle
Calculate and plot spectral energy distribution of main-sequence stars with additional option of blackbody/greybody calculation from mid-IR excess
The code is still under development for general use.

1) Take a look at the Twinkle.ipynb Jupyter notebook to see how to get started.
2) Use sed_paramfile.json to change input parameters
3) Use sample_stardata.txt and ./Docs/sedDatFileDescription.xlsx to set up stellar data and look up definitions and examples of each column used.
4) twinkle.py is the main class called and "inherits" sed.py.

5) sed.py hosts all logistical tools to help calculate the SED: Filter profile management, magnitude/flux conversion,  photosphere generation, etc.

   More documentation is underway


plot_sed.py is a bit more complicated and is being worked on, but basically runs through a list of stars and calculates their photospheric SED as well as a blackbody fit/calculation to any excess associated with specified bands in the sed_paramfile.json file.



Please message me for any comments on improvement or features anyone would like to see.

NextGen and Kurucz photosphere models are currently the only ones incorporated into code. 

References:

1. EMamajek_MSColors.txt from http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
2. NextGen models from Hauschildt, P. H., Allard, F., & Baron, E. 1999, ApJ, 512, 377
3. ATLAS9 modesls from Kurucz, R. L. 1993, yCat, 6039, 0

