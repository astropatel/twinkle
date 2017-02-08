# SED
Calculate and plot spectral energy distribution of main-sequence stars with additional option of blackbody/greybody calculation from mid-IR excess
The code is still under development for general use.

Take a look at the Twinkle.ipynb Jupyter notebook to see how to get started.

plot_sed.py is a bit more complicated and is being worked on, but basically runs through a list of stars and calculates their photospheric SED as well as a blackbody fit/calculation to any excess associated with specified bands in the sed_paramfile.json file.

sed.py hosts all logistical tools to help calculate the SED: Filter profile management, magnitude/flux conversion,  photosphere generation, etc.

More updates underway in the near future.

Please message me for any comments on improvement or features anyone would like to see.

NextGen and Kurucz photosphere models are currently the only ones incorporated into code. 

Ipython Notebook for testing is under development.
