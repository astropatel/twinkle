import os,sys

#=======================================================
#   Supported passbands:
#     
#   Tycho: VT, BT
#   Johnson: UJ, BJ, VJ, RJ, IJ
#   Bessel: UB, BB, VB, RB, IB
#   2MASS: J2M, H2M, Ks2M
#   Spiter/MIPS: MIPS24, MIPS70, MIPS160
#   WISE: W1,W2,W3,W4
#   IRAS: 60, 100
#   Herschel/PACS: 70,100,160
#=======================================================

#////////////////////////////////////////////////////////////////////////////////////
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~FILES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#////////////////////////////////////////////////////////////////////////////////////

#=======================================================
#   DIRECTORY PATH THAT HAS STELLAR INFORMATION
#=======================================================
#file_dir = os.path.join(os.getcwd(),'DebrisDisks','XMatch','FullSky')
file_dir = os.path.join(os.getcwd(),'ExcessStars_07272014_mag2resSNR_fittedCDF_33trim')
#file_dir = os.path.join(os.getcwd(),'WISEExcess_120pcgalplane_thesis')
#file_dir = os.getcwd()

#=======================================================
#    FILE NAME WITH ALL OF THE STELLAR INFORMATION
#=======================================================
#f_read = 'longwave_stars_GPI.txt'
fstar = 'NewWtdExcess_ExcessStars_07272014.txt'
#f_read =   'stars_toplotsed.txt'
#f_read = 'Excesses_120pc_goodstars.txt'

#fstar = 'Excesses_120pc_AllStars_detected.txt'
#f_read = 'jistar.txt'
f_read = os.path.join(file_dir,fstar)

#=======================================================
# FILE NAME AND PATH OF FILE WITH PRE-SELECTED STAR NAMES
#=======================================================
file_select = os.path.join(os.getcwd(),'DebrisDisks','XMatch','FullSky','select_stars.txt')

#=======================================================
#   FILE WITH OPTICAL/NIR EMPIRICAL COLOR RELATIONS
#=======================================================
fileBVStandards = os.path.join(os.getcwd(),'Interpolation_Files','Stellar_colors','EMamajek_MSColors.txt')

#=======================================================
# INITIAL GUESS FOR TEMP AND 10*LOG(g) vs. SPT
#=======================================================
finti_Tg = os.path.join(os.getcwd(),'Interpolation_Files','SED_Init_STG.txt')

#////////////////////////////////////////////////////////////////////////////////////
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    VARIABLES   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#////////////////////////////////////////////////////////////////////////////////////


#=======================================================
#           PLOT PARAMETERS
#=======================================================
# Booleans-- for plotting only one of the two can be on -- cannot plot single and grid at same time
# plot_any overrides plotting -- must be set to True for plotting to occur
plot_any = True
plot_grid = False
plot_single = True
write2file = False
write_SED = False
#redundant
#PLOTTING STUFF
pROW = 4
pCOL= 3
plotsize = (8.78, 10.65)
ylim_up, ylim_low = 10 ** -6, 10 ** -13.8  # PLOT LIMITS IN erg/s/cm^2
xmax = 150
#=======================================================
#          SPECTRUM SAMPLING
#=======================================================
wave_min = 0.2 #microns
wave_max = 4000.0 #microns
gridpts = 1000 #RESOLUTION OF SED

#=======================================================
#  ADOPTION OF W2 AND W3 INTO PHOT MODEL FLAGS
#=======================================================
W3Adapt = True
W2Adapt = True

#=======================================================
#      LAM>30 MICRON ACTIVATION BOOLEAN
#=======================================================
Longwave_Bool = False

#=======================================================
#      PHOTOMETRIC NAME ARRAYS
#=======================================================
# ALL THE PHOTOMETRY THAT IS USED AT ANY POINT IN THE CODE
#mags2use0_original   = ['BJ','VJ','J2M','H2M','Ks2M','W1','W2','W3','W4']
mags2use0_original   = ['BJ','VJ','J2M','H2M','Ks2M','W1','W2','W3','W4']

# PHOTOMETRY THAT IS USED FOR THE STELLAR PORTION OF THE CODE -- PHOTOSPHERIC FIT
mags4Phot0_original  = ['BJ','VJ','J2M','H2M','Ks2M','W1','W2']
#mags4Phot0_original  = ['J2M','H2M','Ks2M']

# PHOTOMETRY THAT IS USED TO SCALE MODELS TO PHOTOMETRIC MEASUREMENTS -- WEIGHTED AVERAGE
mags4scale0_original = ['BJ','VJ','J2M','H2M','Ks2M','W1']
#mags4scale0_original = ['J2M','H2M','Ks2M']

# WHICH PHOTOMETRY IS TO BE USED TO CALCUALTE THE DUST
mags4Dust0  = ['W3','W4']

# BANDS USED TO SCALE SED
scaleSEDbands = ['W1','W2']
#scaleSEDbands = ['W3']

# PHOTOMETRY THAT NEEDS TO BE DISCARDED THAT MIGHT BE SUBJECT TO VARIABILITY -- TYPICALLY LATE K AND M STARS
# CURRENTLY USES COLUMN "NOOPTICAL" FOR EACH STAR TO DETERMINE WHICH STAR THIS APPLIES TO
Remove_RedStars = ['BJ','VJ']


# WHICH FUNCTION TO USE TO MODEL DUST?
Exfunc = 'blackbody'
#exfunc = 'modifiedBB'                    

# WHICH BAND WILL BE USED AS MODIFIED BLACKBODIES PEAK FLUX WAVELENGHT?
modifiedBB_lam0 = 'W3'

mags4DustWFC = ['W3','W4']#,'MIPS24','MIPS70','HPACS70','HPACS100','HPACS160']
#mags4DustWFC = ['W3','W4','MIPS24','MIPS70','HPACS70','HPACS100','HPACS160']

# 99.5 SIGMA CUT OFF
W13_cut = 2.0
W23_cut = 2.0
W12_cut = 4
