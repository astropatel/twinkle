
You can use these filters for the moment. Their corresponding response 
functions are in ./SupportFiles/RSR/.

To specify their use, they need to be added as a string in the "phot" 
section of the sed_paramfile.json

========================================================================
Table of Contents:
------------------
I. File Info
II. Available Filters
III. AB Mag Conversion
IV. Band-relevant info

========================================================================

I. File Info

ALL RSR'S ARE STRICTLY THAT.. ALL LAMBDA*RSR has been divided out

Header portion in each RSR File: #!N0[0,f]/ N1[1,f]/ N2[2,f]/ N3[3,f]/ N4[4,f]/ N5[5,f]/

N0 = Isophotal wavelength (Angstroms)
N1 = vega zero point flux in erg/s/cm^2/Angstrom
N2 = uncertainty in vega zero point flux in erg/s/cm^2/Angstrom
N3 = Isophotal frequency (Hz)
N4 = vega zero point flux in erg/s/cm^2/Hz
N5 = uncertainty in vega zero point flux in erg/s/cm^2/Hz

======================================================================== 

II. Available Filters


    Band Identifiers         Description
    ---------------          ---------------------------

      W1                     WISE BAND 1 (3.4 um)
      W2                     WISE BAND 2 (4.6 um)
      W3                     WISE BAND 3 (~12 um)
      W4                     WISE BAND 4 (~22 um)

      MIPS24                 Spizter/MIPS 24 um 
      MIPS70                 Spizter/MIPS 70 um 
      MIPS160                Spizter/MIPS 160 um

      J2M                    2MASS J
      H2M                    2MASS H
      Ks2M                   2MASS Ks

      UB                     Bessel U
      BB                     Bessel B
      VB                     Bessel V
      RB                     Bessel R
      IB                     Bessel I

      UJ                     Johnson U
      BJ                     Johnson B
      VJ                     Johnson V
      RJ                     Johnson R
      IJ                     Johnson I

      VT                     Tycho V
      BT                     Tycho B
      Hp                     Hipparcos

      IRAS60                 IRAS 60
      IRAS100                IRAS 100
      IRAS25                 IRAS 25
      IRAS12                 IRAS 12

      HPACS70                Herschel/PACS 70
      HPACS100               Herschel/PACS 100
      HPACS160               Herschel/PACS 160

      Akari90                Akari 90 um
      Akari9                 Akari 9 um IRCW
      Akari18                Akari 18 um IRCLW
      
      HNIRC2                 Keck NIRC2 H
      KpNIRC2                Keck NIRC2 Kp
      LpNIRC2                Keck NIRC2 Lp

      MSXA                   MSX A
      MSXC                   MSX C
      MSXD                   MSX D
      MSXE                   MSX E

      IDENIS                 DENIS I BAND
      JDENIS                 DENIS J BAND
      KSDENIS                DENIS Ks BAND

      GGAIA                  GAIA G BAND


========================================================================

III. AB Mag Conversion

AB mag conversion:

http://www.astro.utoronto.ca/~patton/astro/mags.html

http://www.aerith.net/astro/color_conversion.html

http://casa.colorado.edu/~ginsbura/filtersets.htm__

========================================================================

IV. Band Relevant Info

2MASS:

x: wavelength (microns)
y: relative spectral response

U,B,V,I:

x: wavelength (angrstoms)
y: Relative spectral response

Johnson zero points from http://www.stsci.edu/hst/observatory/documents/isrs/scs8.rev.pdf

Taken from
http://spiff.rit.edu/classes/phys440/lectures/filters/filters.html

WISE:
w1,w2,w3,w4
x:wavelength (microns)
y:rsr (not normalized)
# Relative Spectral Response per erg (Equal-Energy); Jarrett et al. 2011
taken from http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html#WISEZMA

The following RSRs are modified RSRs from Brown+2014
The monochromatic flux densities are changed based on
the new vega AB mag = 6.66

    truncated W4:
    Identical to lab RSR except it has zero transmission lambda<21.1 microns

    stretched W4: 
    Identiacal to lab except all the wavelengths have been revised upward by 3.3% (i.e., ∆λ = 0.033λ)
    
    


IRAS:
60
x:wavelength(Angstrom)
y:rsr

SPITZER: 
MIPS: FROM http://irsa.ipac.caltech.edu/data/SPITZER/docs/mips/calibrationfiles/spectralresponse/
Lambda in microns
wavelength is the one gien in MIPS handbook


HERSCHEL:
from SVO profile service
pivot wavelength used
angstroms

AKARI:
from SVO profile service
wavelength in angstrom
Using lam_eff for wavelength: 90
using lam_pivot for wavelength: 18,9


NIRC2:
http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=Keck/NIRC2.H&&mode=browse&gname=Keck&gname2=NIRC2#filter

filter data from NIRC2 website. 
central wavelength and zero point is from SVO 

DENIS, GAIA, MSX:
from SVO profile service
pivot wavelength used.

