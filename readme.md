
<div align="center">
<img src="./docs/source/_static/Logo/twinkle_logo_light.png" alt="Twinkle Logo" width="35%">
</div>

## ReadtheDocs

**Check out the [Twinkle ReadTheDocs Manual](https://twinkle.readthedocs.io/en/latest/)**

## Introduction

**Twinkle** is a Python-based tool that calculates the spectral energy distribution ([SED](https://coolwiki.ipac.caltech.edu/index.php/SED_plots_introduction)) of a stellar source using empirical photometric data and stellar model grids.

*Twinkle* was originally created to help calculate the excess infrared (IR) flux from a star. The presence of an IR excess indicates dust orbiting the star. This dust likely results from the grinding and collisions of asteroids, influenced by a larger planetary object—pointing to the potential for finding planets. You can check out the published papers from my thesis using this code in [Patel, Metchev, and Heinze, 2014](https://iopscience.iop.org/article/10.1088/0067-0049/212/1/10) and [Patel et al., 2017](https://iopscience.iop.org/article/10.3847/1538-3881/153/2/54).

Interested in learning more about debris disks? Check out my [blog post](http://cosmicdiary.org/geminiplanetimager/2015/03/04/debris-disks-searching-for-dust-to-find-planets/).

This code base helps you quickly calculate the temperature and location of the dust to first order by fitting the assumed blackbody or modified blackbody function to the broadband excess emission.

Feel free to fork and contribute!

---

## Available Features

- Model multiple stellar SEDs at once with an easy-to-use stellar input file.
- Late B to K-spectral type modeling support. Additional spectral type modeling is possible but not tested (*yet*).
- Plotting capabilities (link to plotting page upcoming) of the empirical data and the modeled distribution.
- This [Jupyter Notebook](https://github.com/astropatel/twinkle/blob/master/Twinkle_Tutorial.ipynb) provides a quick-start guide to using **Twinkle**.

> **Note:**  
> The code hasn't been tested in a while, so there may be issues still being worked on. If you want to contribute, fork and submit a pull request.

---

## Installation and GitHub

Check out **Twinkle** on [GitHub](https://github.com/astropatel/twinkle).

To install:

```bash
pip install git+https://github.com/astropatel/twinkle.git
```

--- 

## The Logo
I created the Twinkle logo after obtaining baseline ideas from AI image generators, and then drawing it on Inkscape. The top prong, or point, or (what is it called?) of the star image is meant to look like a stellar SED.

---

## The Name
I prefer my projects to have apt names, but then like to backronym it into something non-sensical and fun. While I haven’t landed on one for TWINKLE, here are suggestions from friends and colleagues

The WIse Non-gaseous disK Locator Extraordinaire (by: Scott Barenfeld)

The Wise Infrared Non-gaseous disK Locator Extraordinaire (by: Calen Henderson)

That’s What I Need to Know about teLescopes from Experts (by: Tiffany Meshkat)

---
References:

1. EMamajek_MSColors.txt from http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
2. NextGen models from Hauschildt, P. H., Allard, F., & Baron, E. 1999, ApJ, 512, 377
3. ATLAS9 modesls from Kurucz, R. L. 1993, yCat, 6039, 0

