"""
 mosaic_tools by Rahul I. Patel (ri.patel272@gmail.com)

 Various tools to help plot or analyze monte-carlo results from
 planet detectability simulation.
 Various definitions concerning sorting and creating bins and much more
 There are modules in here for array maniuplation and read/write tools.
 There are also modules here pertaining to fitting.

"""

__author__ = 'Rahul I. Patel'

import scipy
import sys
import random as rnd
import scipy.linalg.blas
import operator
import numpy as np, math as ma
import matplotlib.pyplot as plt
import scipy.interpolate as intp
# import scipy.optimize as opt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from functools import reduce

try:
    from astropy.io import fits
    from astropy.wcs import WCS
except ImportError:
    print('Astropy not detected/installed')


class PlottingTools:
    def triple_axes_dist(ylog=False, xlog=False, xlabel='x', ylabel='y'):
        """
        Sets up plots with 3 axes -- one in center, one on right and one above center
            plot. The purpose is to have a scatter plot or w/e in the center, and two
            distribution plots of the x and y parameters in the side plots.

        Args:
            axScatter: axis object for center scatter plot
            ylog, xlog: booleans to indicate whether x and y axes of center plot are in
                        log scale (base 10)
            xlabel, ylabel : labels for x and y axes labels

            Return
            """
        # pdb.set_trace()
        axScatter = plt.subplot(111)
        axScatter.set_xlabel('%s' % xlabel, fontsize=25)
        axScatter.set_ylabel('%s' % ylabel, fontsize=25)

        divider = make_axes_locatable(axScatter)
        axHistX = divider.append_axes("top", size=2, pad=0.2, sharex=axScatter)
        axHistY = divider.append_axes("right", size=2, pad=0.2,
                                      sharey=axScatter)
        plt.setp(axHistX.get_xticklabels(), visible=False)
        plt.setp(axHistY.get_yticklabels(), visible=False)

        if xlog:
            axScatter.set_xscale('log')
            axHistX.set_xscale('log', nonposy='clip')
        if ylog:
            axScatter.set_yscale('log')
            axHistY.set_yscale('log', nonposy='clip')
        return axScatter, axHistX, axHistY

    def plot_setup(self, axis, gridon=False, minortickson=True,
                   ticklabel_fontsize=20, majortick_width=2.5,
                   minortick_width=1.2, majortick_size=8,
                   minortick_size=5, axes_linewidth=1.5,
                   ytick_direction='in', xtick_direction='in',
                   yaxis_right=False, ylog=False, xlog=False, bold=False,
                   adjust_plot=True):

        """Changes the boring default matplotlib plotting canvas so that it
        looks nice and neat with thicker borders and larger tick marks as well
        as larger fontsizes for the axis labels. Options exist to include or
        exclude the plot grid and minortick mark labels -- set up as boolean
        variables"""

        if gridon:
            axis.grid()
        if minortickson:
            axis.minorticks_on()
        if yaxis_right:
            axis.yaxis.tick_right()

        for line in axis.yaxis.get_majorticklines():
            line.set_markeredgewidth(majortick_width)
        for line in axis.xaxis.get_majorticklines():
            line.set_markeredgewidth(majortick_width)

        for line in axis.xaxis.get_minorticklines():
            line.set_markeredgewidth(minortick_width)
        for line in axis.yaxis.get_minorticklines():
            line.set_markeredgewidth(minortick_width)

        if xlog:
            axis.set_xscale('log', nonposy='clip')
        if ylog:
            axis.set_yscale('log', nonposy='clip')

        # plt.rc('text', usetex=True)
        if bold:
            plt.rc('font', weight='bold')
        plt.rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
        plt.rcParams['mathtext.fontset'] = 'stixsans'
        axis.tick_params(axis='both', which='major',
                         labelsize=ticklabel_fontsize)
        plt.rc("axes", linewidth=axes_linewidth)
        plt.rcParams['xtick.major.size'] = majortick_size
        plt.rcParams['xtick.minor.size'] = minortick_size
        plt.rcParams['ytick.major.size'] = majortick_size
        plt.rcParams['ytick.minor.size'] = minortick_size

        plt.rcParams['xtick.direction'] = xtick_direction
        plt.rcParams['ytick.direction'] = ytick_direction
        if adjust_plot:
            plt.subplots_adjust(left=0.13, bottom=0.13, top=0.95, right=0.97)

        return

    def simpleaxis1(self, ax):
        """This little tool erases the right and top axis lines"""
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        return

    def simpleaxis2(self, ax):
        """This little tool erases the botom and left axis lines"""
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.get_xaxis().tick_top()
        ax.get_yaxis().tick_right()
        return

    def zeroaxes(self, ax):
        ax.spines['left'].set_position('zero')
        ax.spines['right'].set_color('none')
        ax.spines['bottom'].set_position('zero')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_smart_bounds(True)
        ax.spines['bottom'].set_smart_bounds(True)

        return

    def coloraxes(self, ax, color):
        ax.spines['bottom'].set_color('%s' % color)
        ax.spines['top'].set_color('%s' % color)
        ax.spines['right'].set_color('%s' % color)
        ax.spines['left'].set_color('%s' % color)
        return


class ArrayTools:
    """
    A few functions to work with intersection and union of
    numpy arrays
    """

    def intersect_arrays(self, arrays):
        """
        Find intersection of all the arrays in "arrays"
        Returns the sorted, unique values that are in both of the input
        arrays
        """

        N = len(arrays)
        aux = np.array([])

        if N > 1:
            for i in range(N - 1):
                if aux.size == 0:

                    aux = np.intersect1d(arrays[i], arrays[i + 1])
                # print aux
                #                    pdb.set_trace()
                else:

                    aux = np.intersect1d(aux, arrays[i + 1])
                    #                   print aux
                    #                    pdb.set_trace()

        elif N == 1:
            aux = arrays[0]
        else:
            print('No arrays to intersect. Try again.')
            sys.exit()
        return aux

    def union_arrays(self, arrays):

        """
        Find UNION of all the arrays in "arrays"
        Returns the unique, sorted arrays of values that are in either
        of the two input arrays
        """
        N = len(arrays)
        aux = np.array([])

        if N > 1:
            for i in range(N - 1):
                if aux.size == 0:
                    aux = np.union1d(arrays[i], arrays[i + 1])
                else:
                    aux = np.union1d(aux, arrays[i + 1])
        elif N == 1:
            aux = arrays[0]
        else:
            print('No arrays to unionize. Try again.')
            sys.exit()

        return aux

    def dict2list(self, Dict, keys, extra=""):
        """to turn dict into a list.
            Dictionary values which will be converted to numpy array in order
            found in keys.
            dict: Dictionary to be converted to list
            keys: List or array with string values pointing to keys in dict
            Length of keys need not equal length of dict, but len(keys)<= len(dict)
        """
        Dict = Dict
        arr = np.array([])

        for i in range(len(keys)):
            arr = np.append(arr, Dict[keys[i] + extra])

        return arr

    def dictExtract(self, Dict, keys, keySuffix="", newkeySuffix=""):
        """To extract certain values given the input keys
            from the input dictionary and returns a new dictionary
            dict:Dictionary to be sampled from for new dictionary
            keys: keys associated with the input dictionary.
            keySuffix: suffix to be added to each key to access it
            newkeySuffix: if a new key suffix is to be added. Otherwise
            regular keys will be used
            """

        Dict = Dict
        newDict = {}
        for ky in keys:
            newDict[ky + newkeySuffix] = Dict[ky + keySuffix]

        return newDict


class ReadWrite_Tools:
    """
    Read and Write tools once stuff is read

    """

    def create_datadict(self, hnames, data):
        """
        Create a dictionary from a 2D array of data and corresponding header names.

        Args:
            hnames (list): List of header names to use as keys.
            data (array): 2D array where each column corresponds to a header.

        Returns:
            dict: Dictionary with header names as keys and data columns as values.
        """

        p = []
        for j, key in enumerate(hnames):
            p.append((key.strip(), data[j]))
            dat_dict = dict(p)

        return dat_dict

    def create_header(self, list0, more=None, nowrite=None, delimiter='\t '):
        """
        Create a header string from a list of strings, excluding specified items.

        Args:
            list0 (list): List of strings to include in the header.
            more (list, optional): Additional strings to append to the header.
            nowrite (list, optional): Items to exclude from the header.
            delimiter (str, optional): Delimiter used to separate header items.
                Defaults to '\\t'.

        Returns:
            str: Formatted header string with specified items and delimiter.
        """

        listnew = np.array(list0).copy()
        header = ''
        if nowrite is not None:
            nowrite = np.array(nowrite)
        for m in range(len(listnew)):
            if nowrite is None:
                header += listnew[m] + delimiter
            else:

                ind_header = np.where(listnew[m] == nowrite)[0]
                if len(ind_header) == 0:
                    header += listnew[m] + delimiter
                else:
                    pass

        if more is not None:
            for item in more:
                header += item + delimiter
        else:
            pass

        header = header.strip(delimiter)
        header += '\n'
        return header

    def sort_duplicates(self, file, dupcol='object_u', duplicates=None):
        """
        Sort duplicate values in a data file and save a new version with duplicates at the top.

        Args:
            file (str): Path to the input file.
            dupcol (str, optional): Column name used to identify duplicates. Defaults to 'object_u'.
            duplicates (list, optional): Known duplicates. If None, duplicates are identified automatically.

        Returns:
            str: Filename of the new sorted file.
        """

        import collections
        import os
        try:
            from readcol import readcol
        except ImportError:
            print('Readcol not detected/installed')

        # todo: remove readcol
        names, data = readcol(file, names=True)
        # create dictionary with index and column names as val and keys
        # and vice versa
        name_dict2 = dict(enumerate(names))
        name_dict3 = dict(
            list(zip(list(name_dict2.values()), list(name_dict2.keys()))))
        # gather
        dupcolname = data[:, name_dict3[dupcol]]

        # IN CASE NO LIST OF DUPLICATES ARE GIVEN, IT FINDS ALL DUPLICATES AND
        # STORES THOSE NAMES
        if duplicates is None:

            scount = collections.Counter(dupcolname)
            ky_set1, val_set1 = np.array(list(scount.keys())), np.array(
                list(scount.values()))
            ind_dup = np.where(val_set1 > 1)[0]
            duplicate_list = ky_set1[ind_dup]
        else:
            duplicate_list = duplicates

        ind_dup_infile = np.array([])
        for i in range(len(duplicate_list)):
            ind_dup_infilei = np.where(duplicate_list[i] == dupcolname)[0]
            ind_dup_infile = np.append(ind_dup_infile, ind_dup_infilei)

        ind_dup_infile = ind_dup_infile.astype(int)

        # store duplicate data
        dupdata_select = data[ind_dup_infile]

        # Delete duplicate data from original data file
        data = np.delete(data, np.s_[ind_dup_infile], axis=0)

        # Append the sorted duplicate data to top of list
        datanew = np.append(dupdata_select, data, axis=0)

        filebase, file_ext = os.path.splitext(file)

        file2 = filebase + '_2' + file_ext

        names = np.array([names])
        datanew = np.append(names, datanew, axis=0)

        np.savetxt(file2, datanew, fmt='%s', delimiter='\t\t')

        return file2


class FittingTools:
    """
    A collection of utility functions for data fitting tasks.
    """

    def deviates_from_model(self, p0=None, fjac=True, x=None, y=None, err=None,
                            func=None, logx=None, logy=None, loglog=None,
                            **kwargs):
        """
        Returns the deviations (residuals) calculated from an input model function.
        This method is intended to be used by the Levenberg-Marquardt technique in the
        "mpfit.py" module, originally written in IDL by Mark Rivers and Sergey
        Koposov.

        Args:
            p0 (list): Initial parameters for the model to be fit.
            fjac (bool): Flag for partial derivative calculation. Refer to
                MPFIT.py for details.
            x (numpy.ndarray): Input x data (independent variable).
            y (numpy.ndarray): Input y data (dependent variable).
            err (numpy.ndarray, optional): Uncertainty in the y data.
            func (callable): The model function that calculates y from x and
                parameters.
            logx (bool, optional): Flag for applying a logarithmic
                transformation to x.
            logy (bool, optional): Flag for applying a logarithmic
                transformation to y.
            loglog (bool, optional): Flag for applying a logarithmic
                transformation to both x and y.
            **kwargs: Additional parameters to pass to the model function.

        Returns:
            list: A list containing:
                - status (int): Status of the fit, used by the MPFIT.py module.
                - residuals (numpy.ndarray): The residuals, either weighted or unweighted depending on the error input.
        """
        kwargs = kwargs

        # PARAMETERS CAN BE EITHER EXPLICITLY STATED OR WITHIN KWARGS.
        # IF USING MPFIT, PARAMETERS CALLED VIA PARINFO ARE IN KWARGS.
        # OTHERWISE, THEY WILL BE EXPLICIT. EITHER WAY, THEY NEED TO BE EXPLICITLY
        # PASSED TO THE FUNCTION.
        if p0 is not None:
            model = func(x, p0, **kwargs)
        else:
            try:
                p0 = kwargs['p0']
                model = func(x, p0, **kwargs)
                print("Make sure you haven't called parameters twice.")
            except KeyError:
                raise ValueError("No parameters were detected. Try again.")

        status = 0  # needed by mpfit

        #
        if err is not None:
            return ([status, (y - model) / err])
        else:
            return ([status, y - model])

    def poly_nfit(self, x, p):
        """
        Determines the sampled values for a polynomial function whose order is
            determined by the length of the input parameter array.

        Args:
            x (numpy.ndarray): Vector of sampling points.
            p (numpy.ndarray): Vector of parameters for the polynomial (a0, a1, ..., an).

        Returns:
            numpy.ndarray: Sampled values of the polynomial function.
        """
        x, p = np.asarray(x), np.asarray(p)
        x0, p0 = x.copy(), p.copy()

        for i, pi in enumerate(p0):
            i = int(i)
            # MAKE ARRAY FOR FIRST ORDER
            if i == 0:
                y = np.zeros(len(x0)) + pi
            else:
                y += pi * (x0) ** i

        return y

    def get_InitParams(self, lenp):
        """
        Generates an initial set of starting parameters for a fit, depending on
            the number of free parameters in the model.

        Args:
            lenp (int): The number of free parameters in the model.

        Returns:
            list: A list of initial parameters, each set to 0.1.
        """

        return [0.1] * lenp

    def print_poly(self, polyn, xlabel='x', ylabel='y(x)', numformat='%.3f',
                   coeff=None):
        """
        Returns a formatted string for a polynomial function of order N,
        including x and y labels. Optionally, the coefficients can be used to populate the polynomial's parameters.

        Args:
            polyn (int): The order of the polynomial.
            xlabel (str, optional): Label for the x-axis. Defaults to 'x'.
            ylabel (str, optional): Label for the y-axis. Defaults to 'y(x)'.
            numformat (str, optional): String format for displaying the
                coefficients. Defaults to '%.3f'.
            coeff (list, optional): Coefficients for the polynomial.

        Returns:
            str: A string representation of the polynomial.
        """

        sig = {-1: '-', 1: '+'}
        eqstr = r'$%s=' % ylabel

        for i in range(polyn):
            if i == 0:
                eqstr += '%s' + numformat
            elif i == 1:
                eqstr += '%s' + numformat + xlabel
            else:
                eqstr += '%s' + numformat + xlabel + '^%i' % i

        eqstr += '$'
        if coeff is not None:
            eqParams = []
            for p in coeff:
                eqParams.append(sig[np.sign(p)])
                eqParams.append(abs(p))

            eqstr = eqstr % tuple(eqParams)

        else:
            pass

        return eqstr

    def Gauss2d_circle(self, x, p0=None):
        """
        Fits a circular 2D Gaussian centered at (0, 0).

        Args:
            x (numpy.ndarray): The 2D data points (x, y).
            p0 (list, optional): Parameters for the Gaussian model, including
                amplitude and sigma.

        Returns:
            numpy.ndarray: The fitted Gaussian values.
        """
        A, sigma = p0[0], p0[1]
        x0, y0 = x  # this assumes a 2d array
        y = A * np.exp(-(x0 ** 2 + y0 ** 2) / (2 * sigma ** 2))

        return y

    def twoD_Gaussian(self, x, y, amplitude, xo, yo, sigma_x, sigma_y, theta,
                      offset):
        """
        2D Gaussian function for modeling data, adapted from StackOverflow users ali_m and Kokomoking.

        Args:
            x (numpy.ndarray): The x-coordinate data.
            y (numpy.ndarray): The y-coordinate data.
            amplitude (float): The amplitude of the Gaussian.
            xo (float): The x-coordinate of the center.
            yo (float): The y-coordinate of the center.
            sigma_x (float): The standard deviation along the x-axis.
            sigma_y (float): The standard deviation along the y-axis.
            theta (float): The rotation angle of the Gaussian.
            offset (float): The offset value.

        Returns:
            numpy.ndarray: The 2D Gaussian values at the specified coordinates.
        """

        xo = float(xo)
        yo = float(yo)
        a = (np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (
                2 * sigma_y ** 2)
        b = -(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (
                4 * sigma_y ** 2)
        c = (np.sin(theta) ** 2) / (2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (
                2 * sigma_y ** 2)
        g = offset + amplitude * np.exp(
            - (a * ((x - xo) ** 2) + 2 * b * (x - xo) * (y - yo)
               + c * ((y - yo) ** 2)))

        return g.ravel()

    def Gauss_fit(self, x, p0=None):
        """
        Fits a Gaussian function to the data.

        Args:
            x (numpy.ndarray): The input data.
            p0 (list, optional): The parameters for the Gaussian (amplitude, mean, sigma).

        Returns:
            numpy.ndarray: The fitted Gaussian values.
        """

        A, mu, sigma = p0[0], p0[1], p0[2]
        x0 = x
        # y = (A/(sigma*ma.sqrt(2*ma.pi)))*np.exp(-1/2.*((x-mu)/sigma)**2)
        y = A * np.exp(-1 / 2. * ((x - mu) / sigma) ** 2)

        return y

    def exp_fit(self, x, A, b, c):
        """
         Fits an exponential function to the data.

         Args:
            x (numpy.ndarray): The input data.
            A (float): The amplitude.
            b (float): The rate constant.
            c (float): The offset.

         Returns:
            numpy.ndarray: The fitted exponential values.
         """
        x0 = x
        return A * np.exp(b * x0) + c

        return y

    def erf_fit(self, x, p0=None):
        """
        Fits an error function to the data.

        Args:
            x (numpy.ndarray): The input data.
            p0 (list, optional): Parameters for the error function
                (amplitude, mean, sigma).

        Returns:
            numpy.ndarray: The fitted error function values.
        """

        x = np.array(x)
        x0 = x.copy()
        A, mu, sigma = p0.copy()
        y = 0.5 * A * (1. + scipy.special.erf(
            (ma.sqrt(2) / 2) * (x0 - mu) / sigma))
        return y

    def resample_spectrum(self, dataSet1, dataSet2, resample2lowR=False):
        """
        Resamples two input spectra to the same wavelength scale over their common range.

        Args:
            dataSet1, dataSet2 (tuple): Each containing (wavelength, flux) arrays.
            resample2lowR (bool, optional): If True, resamples the high-resolution spectrum to the lower one.

        Returns:
            tuple: A tuple of resampled spectra.
        """
        lam1, flux1 = dataSet1
        lam2, flux2 = dataSet2

        # pdb.set_trace()
        lamMin, lamMax = max(min(lam1), min(lam2)), min(max(lam1), max(lam2))
        index1 = np.where((lam1 >= lamMin) & (lam1 <= lamMax))[0]
        index2 = np.where((lam2 >= lamMin) & (lam2 <= lamMax))[0]
        # index1 = np.where( (lam1>lamMin) & (lam1<lamMax) )[0]
        # index2 = np.where( (lam2>lamMin) & (lam2<lamMax) )[0]

        count1, count2 = len(index1), len(index2)
        try:
            out_lam1 = lam1[index1]
            out_lam2 = lam2[index2]
            out_flx1 = flux1[index1]
            out_flx2 = flux2[index2]
        except IndexError:
            print()
            'Spectra do not overlap in resample_spectrum'
        # if lam1.min()>25349:
        #
        if count1 < count2:
            if not resample2lowR:
                ind = np.where((out_lam2 >= out_lam1.min()) & (
                        out_lam2 <= out_lam1.max()))[0]
                out_lam2 = out_lam2[ind]
                out_flx2 = out_flx2[ind]

                ipolate = intp.interp1d(out_lam1, out_flx1)
                int_flx1 = ipolate(out_lam2)
                out_flx1 = int_flx1
                out_lam1 = out_lam2
            else:  # THIS IS LESS RELIABLE AND UNTESTED -- QUADRATIC PART
                ind = np.where((out_lam1 >= out_lam2.min()) & (
                        out_lam1 <= out_lam2.max()))[0]
                out_lam1 = out_lam1[ind]
                out_flx1 = out_flx1[ind]
                ipolate = intp.interp1d(out_lam2, out_flx2, kind='quadratic')
                int_flx2 = ipolate(out_lam1)
                out_flx2 = int_flx2
                out_lam2 = out_lam1
        else:
            if not resample2lowR:

                ind = np.where((out_lam1 >= out_lam2.min()) & (
                        out_lam1 <= out_lam2.max()))[0]
                out_lam1 = out_lam1[ind]
                out_flx1 = out_flx1[ind]
                ipolate = intp.interp1d(out_lam2, out_flx2)
                int_flx2 = ipolate(out_lam1)
                out_flx2 = int_flx2
                out_lam2 = out_lam1
            else:  # THIS IS LESS RELIABLE AND UNTESTED -- QUADRATIC PART
                ind = np.where((out_lam2 >= out_lam1.min()) & (
                        out_lam2 <= out_lam1.max()))[0]
                out_lam2 = out_lam2[ind]
                out_flx2 = out_flx2[ind]
                ipolate = intp.interp1d(out_lam1, out_flx1, kind='quadratic')
                int_flx1 = ipolate(out_lam2)
                out_flx1 = int_flx1
                out_lam1 = out_lam2

        return ((out_lam1, out_flx1), (out_lam2, out_flx2))

    def resample_model(self, lam, flx, start, end, maxdelta=100.0, pband=None):
        """
        This will resample the input model spectrum to the specified resolution
        between the wavelengths input. If a filter is given, information from the filter
        will be used to supplement the resampling.

        Args:
            lam (arr): array of the wavelength (x-vals) for the spectrum to be resampled
            flx (arr): array of the flux (y-vals) for the spectrum to be resampled
            start,end (float): wavelength bounds for which the resampling should be conducted
            maxdelta (float): maximum difference between wavelegnths tolerated
            pband (obj): passband object

        Returns:
            resampled spectra

        """

        # pdb.set_trace()
        newlam, newflx = [], []

        if len(np.shape(lam)) > 1:
            newlam, newflx = [], []
            for (sublam, subflx) in zip(lam, flx):
                # search for maximum leftmost position in model grid between where filter profile begins
                # and where the largest jump in resolution of model grid
                ind1 = max(np.searchsorted(sublam, start),
                           np.searchsorted(np.diff(sublam), maxdelta))
                # search for position of end of filter profile in model grid
                ind2 = np.searchsorted(sublam, end) + 1
                if ind1 > ind2:
                    newlam.append(sublam)
                    newflx.append(subflx)
                else:
                    # Delta is resolution from ind1 to ind2
                    # only takes first diff calculated -- need to know about other?
                    delta = np.diff(sublam[ind1:ind2])[0]
                    # Split to select array between in1:ind2
                    lams, lamm, laml = np.split(sublam, [ind1, ind2])
                    flxs, flxm, flxl = np.split(subflx, [ind1, ind2])
                    model_interp = intp.interp1d(np.log10(sublam),
                                                 np.log10(subflx))

                    if pband is not None:
                        indpb = np.where((pband.wavelength <= sublam[ind2]) & (
                                pband.wavelength >= sublam[ind1]))[0]
                        lamm = np.unique(
                            np.append(lamm, pband.wavelength[indpb]))
                        # 10 angstrom resolution
                        lamm = np.linspace(lamm[0], lamm[-1],
                                           int(abs(lamm[-1] - lamm[0]) / 10.))

                    else:
                        lamm = np.linspace(lamm[0], lamm[-1],
                                           int(abs(lamm[-1] - lamm[0]) / 10.))

                    new_model_flux = 10 ** model_interp(np.log10(lamm))
                    newsublam = reduce(np.append, [lams, lamm, laml])
                    newsubflx = reduce(np.append, [flxs, new_model_flux, flxl])
                    newlam.append(newsublam)
                    newflx.append(newsubflx)

            maxlenArray = max(list(map(len, newlam)))
            newlam2, newflx2 = [], []
            for k in range(len(newlam)):
                lamThis, flxThis = newlam[k], newflx[k]
                difference = abs(maxlenArray - len(lamThis))
                if difference == 0:
                    pass
                else:
                    ins = np.zeros(difference).astype(int)
                    lamThis = np.append(ins, lamThis)
                    flxThis = np.append(ins, flxThis)
                newlam2.append(lamThis)
                newflx2.append(flxThis)

            newlam, newflx = newlam2, newflx2

        else:

            ind1 = max(np.searchsorted(lam, start),
                       np.searchsorted(np.diff(lam), maxdelta))
            ind2 = np.searchsorted(lam, end)
            if ind1 >= ind2:
                newlam.append(lam)
                newflx.append(flx)
            else:

                delta = np.diff(lam[ind1:ind2])[0]
                lams, lamm, laml = np.split(lam, [ind1, ind2])
                flxs, flxm, flxl = np.split(flx, [ind1, ind2])
                model_interp = intp.interp1d(np.log10(lam), np.log10(flx))

                if pband is not None:
                    indpb = np.where((pband.wavelength <= sublam[ind2]) & (
                            pband.wavelength >= sublam[ind1]))[0]
                    lamm = np.unique(np.append(lamm, pband.wavelength[indpb]))
                    # 10 angstrom resolution
                    lamm = np.linspace(lamm[0], lamm[-1],
                                       int(abs(lamm[-1] - lamm[0]) / 10.))

                else:
                    lamm = np.linspace(lamm[0], lamm[-1],
                                       int(abs(lamm[-1] - lamm[0]) / 10.))

                new_model_flux = 10 ** model_interp(np.log10(lamm))

                newlam = reduce(np.append, [lams, lamm, laml])
                newflx = reduce(np.append, [flxs, new_model_flux, flxl])

        return np.array(newlam), np.array(newflx)


class mpfit:
    """
    Perform Levenberg-Marquardt least-squares minimization, based on MINPACK-1.

    AUTHORS
    The original version of this software, called LMFIT, was written in FORTRAN
    as part of the MINPACK-1 package by XXX.

    Craig Markwardt converted the FORTRAN code to IDL.  The information for the
    IDL version is:
    Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
    craigm@lheamail.gsfc.nasa.gov
    UPDATED VERSIONs can be found on my WEB PAGE:
    http://cow.physics.wisc.edu/~craigm/idl/idl.html

    Mark Rivers created this Python version from Craig's IDL version.
    Mark Rivers, University of Chicago
    Building 434A, Argonne National Laboratory
    9700 South Cass Avenue, Argonne, IL 60439
    rivers@cars.uchicago.edu
    Updated versions can be found at http://cars.uchicago.edu/software

    Sergey Koposov converted the Mark's Python version from Numeric to numpy
    Sergey Koposov, University of Cambridge, Institute of Astronomy,
    Madingley road, CB3 0HA, Cambridge, UK
    koposov@ast.cam.ac.uk
    Updated versions can be found at http://code.google.com/p/astrolibpy/source/browse/trunk/

    DESCRIPTION

    MPFIT uses the Levenberg-Marquardt technique to solve the
    least-squares problem.  In its typical use, MPFIT will be used to
    fit a user-supplied function (the "model") to user-supplied data
    points (the "data") by adjusting a set of parameters.  MPFIT is
    based upon MINPACK-1 (LMDIF.F) by More' and collaborators.

    For example, a researcher may think that a set of observed data
    points is best modelled with a Gaussian curve.  A Gaussian curve is
    parameterized by its mean, standard deviation and normalization.
    MPFIT will, within certain constraints, find the set of parameters
    which best fits the data.  The fit is "best" in the least-squares
    sense; that is, the sum of the weighted squared differences between
    the model and data is minimized.

    The Levenberg-Marquardt technique is a particular strategy for
    iteratively searching for the best fit.  This particular
    implementation is drawn from MINPACK-1 (see NETLIB), and is much faster
    and more accurate than the version provided in the Scientific Python package
    in Scientific.Functions.LeastSquares.
    This version allows upper and lower bounding constraints to be placed on each
    parameter, or the parameter can be held fixed.

    The user-supplied Python function should return an array of weighted
    deviations between model and data.  In a typical scientific problem
    the residuals should be weighted so that each deviate has a
    gaussian sigma of 1.0.  If X represents values of the independent
    variable, Y represents a measurement for each value of X, and ERR
    represents the error in the measurements, then the deviates could
    be calculated as follows:

    DEVIATES = (Y - F(X)) / ERR

    where F is the analytical function representing the model.  You are
    recommended to use the convenience functions MPFITFUN and
    MPFITEXPR, which are driver functions that calculate the deviates
    for you.  If ERR are the 1-sigma uncertainties in Y, then

    TOTAL( DEVIATES^2 )

    will be the total chi-squared value.  MPFIT will minimize the
    chi-square value.  The values of X, Y and ERR are passed through
    MPFIT to the user-supplied function via the FUNCTKW keyword.

    Simple constraints can be placed on parameter values by using the
    PARINFO keyword to MPFIT.  See below for a description of this
    keyword.

    MPFIT does not perform more general optimization tasks.  See TNMIN
    instead.  MPFIT is customized, based on MINPACK-1, to the
    least-squares minimization problem.

    USER FUNCTION

    The user must define a function which returns the appropriate
    values as specified above.  The function should return the weighted
    deviations between the model and the data.  It should also return a status
    flag and an optional partial derivative array.  For applications which
    use finite-difference derivatives -- the default -- the user
    function should be declared in the following way:

    def myfunct(p, fjac=None, x=None, y=None, err=None)
    # Parameter values are passed in "p"
    # If fjac==None then partial derivatives should not be
    # computed.  It will always be None if MPFIT is called with default
    # flag.
    model = F(x, p)
    # Non-negative status value means MPFIT should continue, negative means
    # stop the calculation.
    status = 0
    return([status, (y-model)/err]

    See below for applications with analytical derivatives.

    The keyword parameters X, Y, and ERR in the example above are
    suggestive but not required.  Any parameters can be passed to
    MYFUNCT by using the functkw keyword to MPFIT.  Use MPFITFUN and
    MPFITEXPR if you need ideas on how to do that.  The function *must*
    accept a parameter list, P.

    In general there are no restrictions on the number of dimensions in
    X, Y or ERR.  However the deviates *must* be returned in a
    one-dimensional Numeric array of type Float.

    User functions may also indicate a fatal error condition using the
    status return described above. If status is set to a number between
    -15 and -1 then MPFIT will stop the calculation and return to the caller.


                                ANALYTIC DERIVATIVES

    In the search for the best-fit solution, MPFIT by default
    calculates derivatives numerically via a finite difference
    approximation.  The user-supplied function need not calculate the
    derivatives explicitly.  However, if you desire to compute them
    analytically, then the AUTODERIVATIVE=0 keyword must be passed to MPFIT.
    As a practical matter, it is often sufficient and even faster to allow
    MPFIT to calculate the derivatives numerically, and so
    AUTODERIVATIVE=0 is not necessary.

    If AUTODERIVATIVE=0 is used then the user function must check the parameter
    FJAC, and if FJAC!=None then return the partial derivative array in the
    return list.

    where FGRAD(x, p, i) is a user function which must compute the
    derivative of the model with respect to parameter P[i] at X.  When
    finite differencing is used for computing derivatives (ie, when
    AUTODERIVATIVE=1), or when MPFIT needs only the errors but not the
    derivatives the parameter FJAC=None.

    Derivatives should be returned in the PDERIV array. PDERIV should be an m x
    n array, where m is the number of data points and n is the number
    of parameters.  dp[i,j] is the derivative at the ith point with
    respect to the jth parameter.

    The derivatives with respect to fixed parameters are ignored; zero
    is an appropriate value to insert for those derivatives.  Upon
    input to the user function, FJAC is set to a vector with the same
    length as P, with a value of 1 for a parameter which is free, and a
    value of zero for a parameter which is fixed (and hence no
    derivative needs to be calculated).

    If the data is higher than one dimensional, then the *last*
    dimension should be the parameter dimension.  Example: fitting a
    50x50 image, "dp" should be 50x50xNPAR.


    CONSTRAINING PARAMETER VALUES WITH THE PARINFO KEYWORD

    The behavior of MPFIT can be modified with respect to each
    parameter to be fitted.  A parameter value can be fixed; simple
    boundary constraints can be imposed; limitations on the parameter
    changes can be imposed; properties of the automatic derivative can
    be modified; and parameters can be tied to one another.

    These properties are governed by the PARINFO structure, which is
    passed as a keyword parameter to MPFIT.

    PARINFO should be a list of dictionaries, one list entry for each parameter.
    Each parameter is associated with one element of the array, in
    numerical order.  The dictionary can have the following keys
    (none are required, keys are case insensitive):

    'value' - the starting parameter value (but see the START_PARAMS
         parameter for more information).

    'fixed' - a boolean value, whether the parameter is to be held
         fixed or not.  Fixed parameters are not varied by
         MPFIT, but are passed on to MYFUNCT for evaluation.

    'limited' - a two-element boolean array.  If the first/second
           element is set, then the parameter is bounded on the
           lower/upper side.  A parameter can be bounded on both
           sides.  Both LIMITED and LIMITS must be given
           together.

    'limits' - a two-element float array.  Gives the
          parameter limits on the lower and upper sides,
          respectively.  Zero, one or two of these values can be
          set, depending on the values of LIMITED.  Both LIMITED
          and LIMITS must be given together.

    'parname' - a string, giving the name of the parameter.  The
           fitting code of MPFIT does not use this tag in any
           way.  However, the default iterfunct will print the
           parameter name if available.

    'step' - the step size to be used in calculating the numerical
        derivatives.  If set to zero, then the step size is
        computed automatically.  Ignored when AUTODERIVATIVE=0.

    'mpside' - the sidedness of the finite difference when computing
          numerical derivatives.  This field can take four
          values:

             0 - one-sided derivative computed automatically
             1 - one-sided derivative (f(x+h) - f(x)  )/h
            -1 - one-sided derivative (f(x)   - f(x-h))/h
             2 - two-sided derivative (f(x+h) - f(x-h))/(2*h)

         Where H is the STEP parameter described above.  The
         "automatic" one-sided derivative method will chose a
         direction for the finite difference which does not
         violate any constraints.  The other methods do not
         perform this check.  The two-sided method is in
         principle more precise, but requires twice as many
         function evaluations.  Default: 0.

    'mpmaxstep' - the maximum change to be made in the parameter
             value.  During the fitting process, the parameter
             will never be changed by more than this value in
             one iteration.

             A value of 0 indicates no maximum.  Default: 0.

    'tied' - a string expression which "ties" the parameter to other
        free or fixed parameters.  Any expression involving
        constants and the parameter array P are permitted.
        Example: if parameter 2 is always to be twice parameter
        1 then use the following: parinfo(2).tied = '2 * p(1)'.
        Since they are totally constrained, tied parameters are
        considered to be fixed; no errors are computed for them.
        [ NOTE: the PARNAME can't be used in expressions. ]

    'mpprint' - if set to 1, then the default iterfunct will print the
           parameter value.  If set to 0, the parameter value
           will not be printed.  This tag can be used to
           selectively print only a few parameter values out of
           many.  Default: 1 (all parameters printed)


    Future modifications to the PARINFO structure, if any, will involve
    adding dictionary tags beginning with the two letters "MP".
    Therefore programmers are urged to avoid using tags starting with
    the same letters; otherwise they are free to include their own
    fields within the PARINFO structure, and they will be ignored.

    PARINFO Example:
    parinfo = [{'value':0., 'fixed':0, 'limited':[0,0], 'limits':[0.,0.]}
                                                     for i in range(5)]
    parinfo[0]['fixed'] = 1
    parinfo[4]['limited'][0] = 1
    parinfo[4]['limits'][0]  = 50.
    values = [5.7, 2.2, 500., 1.5, 2000.]
    for i in range(5): parinfo[i]['value']=values[i]

    A total of 5 parameters, with starting values of 5.7,
    2.2, 500, 1.5, and 2000 are given.  The first parameter
    is fixed at a value of 5.7, and the last parameter is
    constrained to be above 50.


    EXAMPLE

    import mpfit
    import np.oldnumeric as Numeric
    x = arange(100, float)
    p0 = [5.7, 2.2, 500., 1.5, 2000.]
    y = ( p[0] + p[1]*[x] + p[2]*[x**2] + p[3]*sqrt(x) + p[4]*log(x))
    fa = {'x':x, 'y':y, 'err':err}
    m = mpfit('myfunct', p0, functkw=fa)
    print 'status = ', m.status
    if (m.status <= 0): print 'error message = ', m.errmsg
    print 'parameters = ', m.params

    Minimizes sum of squares of MYFUNCT.  MYFUNCT is called with the X,
    Y, and ERR keyword parameters that are given by FUNCTKW.  The
    results can be obtained from the returned object m.


    THEORY OF OPERATION

    There are many specific strategies for function minimization.  One
    very popular technique is to use function gradient information to
    realize the local structure of the function.  Near a local minimum
    the function value can be taylor expanded about x0 as follows:

    f(x) = f(x0) + f'(x0) . (x-x0) + (1/2) (x-x0) . f''(x0) . (x-x0)
    Order   0th          1st                     2nd

    Here f'(x) is the gradient vector of f at x, and f''(x) is the
    Hessian matrix of second derivatives of f at x.  The vector x is
    the set of function parameters, not the measured data vector.  One
    can find the minimum of f, f(xm) using Newton's method, and
    arrives at the following linear equation:

    f''(x0) . (xm-x0) = - f'(x0)                            (2)

    If an inverse can be found for f''(x0) then one can solve for
    (xm-x0), the step vector from the current position x0 to the new
    projected minimum.  Here the problem has been linearized (ie, the
    gradient information is known to first order).  f''(x0) is
    symmetric n x n matrix, and should be positive definite.

    The Levenberg - Marquardt technique is a variation on this theme.
    It adds an additional diagonal term to the equation which may aid the
    convergence properties:

    (f''(x0) + nu I) . (xm-x0) = -f'(x0)                  (2a)

    where I is the identity matrix.  When nu is large, the overall
    matrix is diagonally dominant, and the iterations follow steepest
    descent.  When nu is small, the iterations are quadratically
    convergent.

    In principle, if f''(x0) and f'(x0) are known then xm-x0 can be
    determined.  However the Hessian matrix is often difficult or
    impossible to compute.  The gradient f'(x0) may be easier to
    compute, if even by finite difference techniques.  So-called
    quasi-Newton techniques attempt to successively estimate f''(x0)
    by building up gradient information as the iterations proceed.

    In the least squares problem there are further simplifications
    which assist in solving eqn (2).  The function to be minimized is
    a sum of squares:

    f = Sum(hi^2)                                         (3)

    where hi is the ith residual out of m residuals as described
    above.  This can be substituted back into eqn (2) after computing
    the derivatives:

    f'  = 2 Sum(hi  hi')
    f'' = 2 Sum(hi' hj') + 2 Sum(hi hi'')                (4)

    If one assumes that the parameters are already close enough to a
    minimum, then one typically finds that the second term in f'' is
    negligible [or, in any case, is too difficult to compute].  Thus,
    equation (2) can be solved, at least approximately, using only
    gradient information.

    In matrix notation, the combination of eqns (2) and (4) becomes:

    hT' . h' . dx = - hT' . h                          (5)

    Where h is the residual vector (length m), hT is its transpose, h'
    is the Jacobian matrix (dimensions n x m), and dx is (xm-x0).  The
    user function supplies the residual vector h, and in some cases h'
    when it is not found by finite differences (see MPFIT_FDJAC2,
    which finds h and hT').  Even if dx is not the best absolute step
    to take, it does provide a good estimate of the best *direction*,
    so often a line minimization will occur along the dx vector
    direction.

    The method of solution employed by MINPACK is to form the Q . R
    factorization of h', where Q is an orthogonal matrix such that QT .
    Q = I, and R is upper right triangular.  Using h' = Q . R and the
    ortogonality of Q, eqn (5) becomes

    (RT . QT) . (Q . R) . dx = - (RT . QT) . h
    RT . R . dx = - RT . QT . h          (6)
    R . dx = - QT . h

    where the last statement follows because R is upper triangular.
    Here, R, QT and h are known so this is a matter of solving for dx.
    The routine MPFIT_QRFAC provides the QR factorization of h, with
    pivoting, and MPFIT_QRSOLV provides the solution for dx.


    REFERENCES

    MINPACK-1, Jorge More', available from netlib (www.netlib.org).
    "Optimization Software Guide," Jorge More' and Stephen Wright, SIAM,
    *Frontiers in Applied Mathematics*, Number 14.
    More', Jorge J., "The Levenberg-Marquardt Algorithm: Implementation and
    Theory," in *Numerical Analysis*, ed. Watson, G. A., Lecture Notes
    in Mathematics 630, Springer-Verlag, 1977.


    MODIFICATION HISTORY

    Translated from MINPACK-1 in FORTRAN, Apr-Jul 1998, CM
    Copyright (C) 1997-2002, Craig Markwardt
    This software is provided as is without any warranty whatsoever.
    Permission to use, copy, modify, and distribute modified or
    unmodified copies is granted, provided this copyright and disclaimer
    are included unchanged.

    Translated from MPFIT (Craig Markwardt's IDL package) to Python,
    August, 2002.  Mark Rivers
    Converted from Numeric to numpy (Sergey Koposov, July 2008)
"""

    blas_enorm32, = scipy.linalg.blas.get_blas_funcs(['nrm2'], np.array([0],
                                                                        dtype=np.float32))
    blas_enorm64, = scipy.linalg.blas.get_blas_funcs(['nrm2'], np.array([0],
                                                                        dtype=np.float64))

    def __init__(self, fcn, xall=None, functkw={}, parinfo=None,
                 ftol=1.e-10, xtol=1.e-10, gtol=1.e-10,
                 damp=0., maxiter=2000, factor=100., nprint=1,
                 iterfunct='default', iterkw={}, nocovar=0,
                 rescale=0, autoderivative=1, quiet=0,
                 diag=None, epsfcn=None, debug=0):
        """
        Initialize the optimizer with user-defined parameters.

        Args:
            fcn (callable): The function to be minimized. The function should return
                the weighted deviations between the model and the data.
            xall (array, optional): An array of starting values for each of the
                parameters of the model. The number of parameters should be fewer
                than the number of measurements. This parameter is optional if the
                `parinfo` keyword is used.
            functkw (dict, optional): A dictionary containing parameters to be
                passed to the user-supplied function specified by `fcn` via the
                keyword dictionary mechanism. Default is {} (no extra parameters).
            parinfo (list of dicts, optional): A mechanism for more sophisticated
                constraints to be placed on parameter values. Default is None (all
                parameters are free and unconstrained).
            ftol (float, optional): A nonnegative variable that measures the
                desired relative error in the sum of squares. Default is 1E-10.
            xtol (float, optional): A nonnegative variable that measures the
                desired relative error in the approximate solution. Default is 1E-10.
            gtol (float, optional): A nonnegative variable that measures the
                orthogonality desired between the function vector and the columns
                of the Jacobian. Default is 1E-10.
            damp (float, optional): A scalar number indicating the cut-off value
                of residuals where damping will occur. Default is 0 (no damping).
            maxiter (int, optional): The maximum number of iterations to perform.
                Default is 2000.
            factor (float, optional): A factor that influences the step size.
                Default is 100.
            nprint (int, optional): The frequency with which `iterfunct` is called.
                Default is 1.
            iterfunct (callable, optional): A function to be called upon each
                `nprint` iteration of the MPFIT routine. Default is 'default'.
            iterkw (dict, optional): Keyword arguments to be passed to `iterfunct`.
                Default is {} (no arguments passed).
            nocovar (int, optional): Set this keyword to prevent the calculation
                of the covariance matrix before returning. Default is 0 (covariance
                matrix is returned).
            rescale (int, optional): Flag for rescaling parameters. Default is 0.
            autoderivative (int, optional): If set, derivatives of the function
                will be computed automatically. Default is 1 (set).
            quiet (int, optional): Set this keyword when no textual output
                should be printed by MPFIT. Default is 0.
            diag (array, optional): An array for parameter scaling. Default is None.
            epsfcn (float, optional): A parameter for finite differencing.
                Default is None.
            debug (int, optional): Debugging flag. Default is 0.

        Returns:
            None: This function initializes the instance variables for the class,
            which will hold the optimization results as attributes.

            .params: The current set of model parameters.
            .niter: The number of iterations completed.
            .covar: The covariance matrix for the set of parameters. It is NxN,
                    where N is the number of parameters. If `nocovar` is set,
                    this will be None.
            .perror: The formal 1-sigma errors in each parameter, computed from
                     the covariance matrix. If a parameter is held fixed or if it
                     touches a boundary, the error is reported as zero.
            .status: An integer status code representing the result of the
                     optimization. Possible values include:
                -16: A parameter or function value has become infinite, usually
                     due to numerical overflow.
                -15 to -1: Error codes that either the user-supplied
                     function or `iterfunct` may return to terminate the fitting
                     process.
                0: Improper input parameters.
                1: Both actual and predicted relative reductions in the sum
                     of squares are at most `ftol`.
                2: Relative error between two consecutive iterates is at most `xtol`.
                3: Conditions for status = 1 and status = 2 both hold.
                4: The cosine of the angle between `fvec` and any column of
                     the Jacobian is at most `gtol` in absolute value.
                5: The maximum number of iterations has been reached.
                6: `ftol` is too small; no further reduction in the sum of squares is
                     possible.
                7: `xtol` is too small; no further improvement in the approximate
                     solution is possible.
                8: `gtol` is too small; `fvec` is orthogonal to the columns
                     of the Jacobian to machine precision.
            .errmsg: A string error or warning message, providing additional context
                     for the status of the optimization.
            .nfev: The number of calls to the user-defined function performed.
            .damp: The damping parameter.

        """
        self.niter = 0
        self.params = None
        self.covar = None
        self.perror = None
        self.status = 0  # Invalid input flag set while we check inputs
        self.debug = debug
        self.errmsg = ''
        self.nfev = 0
        self.damp = damp
        self.dof = 0

        if fcn == None:
            self.errmsg = "Usage: parms = mpfit('myfunt', ... )"
            return

        if iterfunct == 'default':
            iterfunct = self.defiter

        # Parameter damping doesn't work when user is providing their own
        # gradients.
        if (self.damp != 0) and (autoderivative == 0):
            self.errmsg = 'ERROR: keywords DAMP and AUTODERIVATIVE are mutually exclusive'
            return

        # Parameters can either be stored in parinfo, or x. x takes precedence if it exists
        if (xall is None) and (parinfo is None):
            self.errmsg = 'ERROR: must pass parameters in P or PARINFO'
            return

        # Be sure that PARINFO is of the right type
        if parinfo is not None:
            if type(parinfo) != list:
                self.errmsg = 'ERROR: PARINFO must be a list of dictionaries.'
                return
            else:
                if type(parinfo[0]) != dict:
                    self.errmsg = 'ERROR: PARINFO must be a list of dictionaries.'
                    return
            if ((xall is not None) and (len(xall) != len(parinfo))):
                self.errmsg = 'ERROR: number of elements in PARINFO and P must agree'
                return

        # If the parameters were not specified at the command line, then
        # extract them from PARINFO
        if xall is None:
            xall = self.parinfo(parinfo, 'value')
            if xall is None:
                self.errmsg = 'ERROR: either P or PARINFO(*)["value"] must be supplied.'
                return

        # Make sure parameters are numpy arrays
        xall = np.asarray(xall)
        # In the case if the xall is not float or if is float but has less
        # than 64 bits we do convert it into double
        if xall.dtype.kind != 'f' or xall.dtype.itemsize <= 4:
            xall = xall.astype(np.float)

        npar = len(xall)
        self.fnorm = -1.
        fnorm1 = -1.

        # TIED parameters?
        ptied = self.parinfo(parinfo, 'tied', default='', n=npar)
        self.qanytied = 0
        for i in range(npar):
            ptied[i] = ptied[i].strip()
            if ptied[i] != '':
                self.qanytied = 1
        self.ptied = ptied

        # FIXED parameters ?
        pfixed = self.parinfo(parinfo, 'fixed', default=0, n=npar)
        pfixed = (pfixed == 1)
        for i in range(npar):
            pfixed[i] = pfixed[i] or (ptied[
                                          i] != '')  # Tied parameters are also effectively fixed

        # Finite differencing step, absolute and relative, and sidedness of deriv.
        step = self.parinfo(parinfo, 'step', default=0., n=npar)
        dstep = self.parinfo(parinfo, 'relstep', default=0., n=npar)
        dside = self.parinfo(parinfo, 'mpside', default=0, n=npar)

        # Maximum and minimum steps allowed to be taken in one iteration
        maxstep = self.parinfo(parinfo, 'mpmaxstep', default=0., n=npar)
        minstep = self.parinfo(parinfo, 'mpminstep', default=0., n=npar)
        qmin = minstep != 0
        qmin[:] = False  # Remove minstep for now!!
        qmax = maxstep != 0
        if np.any(qmin & qmax & (maxstep < minstep)):
            self.errmsg = 'ERROR: MPMINSTEP is greater than MPMAXSTEP'
            return
        wh = (np.nonzero((qmin != 0.) | (qmax != 0.)))[0]
        qminmax = len(wh > 0)

        # Finish up the free parameters
        ifree = (np.nonzero(pfixed != 1))[0]
        nfree = len(ifree)
        if nfree == 0:
            self.errmsg = 'ERROR: no free parameters'
            return

        # Compose only VARYING parameters
        self.params = xall.copy()  # self.params is the set of parameters to be returned
        x = self.params[ifree]  # x is the set of free parameters

        # LIMITED parameters ?
        limited = self.parinfo(parinfo, 'limited', default=[0, 0], n=npar)
        limits = self.parinfo(parinfo, 'limits', default=[0., 0.], n=npar)
        if (limited is not None) and (limits is not None):
            # Error checking on limits in parinfo
            if np.any((limited[:, 0] & (xall < limits[:, 0])) |
                      (limited[:, 1] & (xall > limits[:, 1]))):
                self.errmsg = 'ERROR: parameters are not within PARINFO limits'
                return
            if np.any((limited[:, 0] & limited[:, 1]) &
                      (limits[:, 0] >= limits[:, 1]) &
                      (pfixed == 0)):
                self.errmsg = 'ERROR: PARINFO parameter limits are not consistent'
                return

            # Transfer structure values to local variables
            qulim = (limited[:, 1])[ifree]
            ulim = (limits[:, 1])[ifree]
            qllim = (limited[:, 0])[ifree]
            llim = (limits[:, 0])[ifree]

            if np.any((qulim != 0.) | (qllim != 0.)):
                qanylim = 1
            else:
                qanylim = 0
        else:
            # Fill in local variables with dummy values
            qulim = np.zeros(nfree)
            ulim = x * 0.
            qllim = qulim
            llim = x * 0.
            qanylim = 0

        n = len(x)
        # Check input parameters for errors
        if (n < 0) or (ftol <= 0) or (xtol <= 0) or (gtol <= 0) \
                or (maxiter < 0) or (factor <= 0):
            self.errmsg = 'ERROR: input keywords are inconsistent'
            return

        if rescale != 0:
            self.errmsg = 'ERROR: DIAG parameter scales are inconsistent'
            if len(diag) < n:
                return
            if np.any(diag <= 0):
                return
            self.errmsg = ''

        [self.status, fvec] = self.call(fcn, self.params, functkw)

        if self.status < 0:
            self.errmsg = 'ERROR: first call to "' + str(fcn) + '" failed'
            return
        # If the returned fvec has more than four bits I assume that we have
        # double precision
        # It is important that the machar is determined by the precision of
        # the returned value, not by the precision of the input array
        if np.array([fvec]).dtype.itemsize > 4:
            self.machar = machar(double=1)
            self.blas_enorm = mpfit.blas_enorm64
        else:
            self.machar = machar(double=0)
            self.blas_enorm = mpfit.blas_enorm32
        machep = self.machar.machep

        m = len(fvec)
        if m < n:
            self.errmsg = 'ERROR: number of parameters must not exceed data'
            return
        self.dof = m - nfree
        self.fnorm = self.enorm(fvec)

        # Initialize Levelberg-Marquardt parameter and iteration counter

        par = 0.
        self.niter = 1
        qtf = x * 0.
        self.status = 0

        # Beginning of the outer loop

        while (1):

            # If requested, call fcn to enable printing of iterates
            self.params[ifree] = x
            if self.qanytied:
                self.params = self.tie(self.params, ptied)

            if (nprint > 0) and (iterfunct is not None):
                if ((self.niter - 1) % nprint) == 0:
                    mperr = 0
                    xnew0 = self.params.copy()

                    dof = np.max([len(fvec) - len(x), 0])
                    status = iterfunct(fcn, self.params, self.niter,
                                       self.fnorm ** 2,
                                       functkw=functkw, parinfo=parinfo,
                                       quiet=quiet,
                                       dof=dof, **iterkw)
                    if status is not None:
                        self.status = status

                    # Check for user termination
                    if self.status < 0:
                        self.errmsg = 'WARNING: premature termination by ' + str(
                            iterfunct)
                        return

                    # If parameters were changed (grrr..) then re-tie
                    if np.max(np.abs(xnew0 - self.params)) > 0:
                        if self.qanytied:
                            self.params = self.tie(self.params, ptied)
                        x = self.params[ifree]

            # Calculate the jacobian matrix
            self.status = 2
            catch_msg = 'calling MPFIT_FDJAC2'
            #
            fjac = self.fdjac2(fcn, x, fvec, step, qulim, ulim, dside,
                               epsfcn=epsfcn,
                               autoderivative=autoderivative, dstep=dstep,
                               functkw=functkw, ifree=ifree, xall=self.params)
            if fjac is None:
                self.errmsg = 'WARNING: premature termination by FDJAC2'
                return

            # Determine if any of the parameters are pegged at the limits
            if qanylim:
                catch_msg = 'zeroing derivatives of pegged parameters'
                whlpeg = (np.nonzero(qllim & (x == llim)))[0]
                nlpeg = len(whlpeg)
                whupeg = (np.nonzero(qulim & (x == ulim)))[0]
                nupeg = len(whupeg)
                # See if any "pegged" values should keep their derivatives
                if nlpeg > 0:
                    # Total derivative of sum wrt lower pegged parameters
                    for i in range(nlpeg):
                        sum0 = sum(fvec * fjac[:, whlpeg[i]])
                        if sum0 > 0:
                            fjac[:, whlpeg[i]] = 0
                if nupeg > 0:
                    # Total derivative of sum wrt upper pegged parameters
                    for i in range(nupeg):
                        sum0 = sum(fvec * fjac[:, whupeg[i]])
                        if sum0 < 0:
                            fjac[:, whupeg[i]] = 0

            # Compute the QR factorization of the jacobian
            [fjac, ipvt, wa1, wa2] = self.qrfac(fjac, pivot=1)

            # On the first iteration if "diag" is unspecified, scale
            # according to the norms of the columns of the initial jacobian
            catch_msg = 'rescaling diagonal elements'
            if self.niter == 1:
                if (rescale == 0) or (len(diag) < n):
                    diag = wa2.copy()
                    diag[diag == 0] = 1.

                # On the first iteration, calculate the norm of the scaled x
                # and initialize the step bound delta
                wa3 = diag * x
                xnorm = self.enorm(wa3)
                delta = factor * xnorm
                if delta == 0.:
                    delta = factor

            # Form (q transpose)*fvec and store the first n components in qtf
            catch_msg = 'forming (q transpose)*fvec'
            wa4 = fvec.copy()
            for j in range(n):
                lj = ipvt[j]
                temp3 = fjac[j, lj]
                if temp3 != 0:
                    fj = fjac[j:, lj]
                    wj = wa4[j:]
                    # *** optimization wa4(j:*)
                    wa4[j:] = wj - fj * sum(fj * wj) / temp3
                fjac[j, lj] = wa1[j]
                qtf[j] = wa4[j]
            # From this point on, only the square matrix, consisting of the
            # triangle of R, is needed.
            fjac = fjac[0:n, 0:n]
            fjac.shape = [n, n]
            temp = fjac.copy()
            for i in range(n):
                temp[:, i] = fjac[:, ipvt[i]]
            fjac = temp.copy()

            # Check for overflow.  This should be a cheap test here since FJAC
            # has been reduced to a (small) square matrix, and the test is
            # O(N^2).
            # wh = where(finite(fjac) EQ 0, ct)
            # if ct GT 0 then goto, FAIL_OVERFLOW

            # Compute the norm of the scaled gradient
            catch_msg = 'computing the scaled gradient'
            gnorm = 0.
            if self.fnorm != 0:
                for j in range(n):
                    l = ipvt[j]
                    if wa2[l] != 0:
                        sum0 = sum(fjac[0:j + 1, j] * qtf[0:j + 1]) / self.fnorm
                        gnorm = np.max([gnorm, np.abs(sum0 / wa2[l])])

            # Test for convergence of the gradient norm
            if gnorm <= gtol:
                self.status = 4
                break
            if maxiter == 0:
                self.status = 5
                break

            # Rescale if necessary
            if rescale == 0:
                diag = np.choose(diag > wa2, (wa2, diag))

            # Beginning of the inner loop
            while (1):

                # Determine the levenberg-marquardt parameter
                catch_msg = 'calculating LM parameter (MPFIT_)'
                [fjac, par, wa1, wa2] = self.lmpar(fjac, ipvt, diag, qtf,
                                                   delta, wa1, wa2, par=par)
                # Store the direction p and x+p. Calculate the norm of p
                wa1 = -wa1

                if (qanylim == 0) and (qminmax == 0):
                    # No parameter limits, so just move to new position WA2
                    alpha = 1.
                    wa2 = x + wa1

                else:

                    # Respect the limits.  If a step were to go out of bounds, then
                    # we should take a step in the same direction but shorter distance.
                    # The step should take us right to the limit in that case.
                    alpha = 1.

                    if qanylim:
                        # Do not allow any steps out of bounds
                        catch_msg = 'checking for a step out of bounds'
                        if nlpeg > 0:
                            wa1[whlpeg] = np.clip(wa1[whlpeg], 0., np.max(wa1))
                        if nupeg > 0:
                            wa1[whupeg] = np.clip(wa1[whupeg], np.min(wa1), 0.)

                        dwa1 = np.abs(wa1) > machep
                        whl = (np.nonzero(
                            ((dwa1 != 0.) & qllim) & ((x + wa1) < llim)))[0]
                        if len(whl) > 0:
                            t = ((llim[whl] - x[whl]) /
                                 wa1[whl])
                            alpha = np.min([alpha, np.min(t)])
                        whu = (np.nonzero(
                            ((dwa1 != 0.) & qulim) & ((x + wa1) > ulim)))[0]
                        if len(whu) > 0:
                            t = ((ulim[whu] - x[whu]) /
                                 wa1[whu])
                            alpha = np.min([alpha, np.min(t)])

                    # Obey any max step values.
                    if qminmax:
                        nwa1 = wa1 * alpha
                        whmax = (np.nonzero((qmax != 0.) & (maxstep > 0)))[0]
                        if len(whmax) > 0:
                            mrat = np.max(np.abs(nwa1[whmax]) /
                                          np.abs(maxstep[ifree[whmax]]))
                            if mrat > 1:
                                alpha = alpha / mrat

                    # Scale the resulting vector
                    wa1 = wa1 * alpha
                    wa2 = x + wa1

                    # Adjust the final output values.  If the step put us exactly
                    # on a boundary, make sure it is exact.
                    sgnu = (ulim >= 0) * 2. - 1.
                    sgnl = (llim >= 0) * 2. - 1.
                    # Handles case of
                    #        ... nonzero *LIM ... ...zero * LIM
                    ulim1 = ulim * (1 - sgnu * machep) - (ulim == 0) * machep
                    llim1 = llim * (1 + sgnl * machep) + (llim == 0) * machep
                    wh = (np.nonzero((qulim != 0) & (wa2 >= ulim1)))[0]
                    if len(wh) > 0:
                        wa2[wh] = ulim[wh]
                    wh = (np.nonzero((qllim != 0.) & (wa2 <= llim1)))[0]
                    if len(wh) > 0:
                        wa2[wh] = llim[wh]
                # endelse
                wa3 = diag * wa1
                pnorm = self.enorm(wa3)

                # On the first iteration, adjust the initial step bound
                if self.niter == 1:
                    delta = np.min([delta, pnorm])

                self.params[ifree] = wa2

                # Evaluate the function at x+p and calculate its norm
                mperr = 0
                catch_msg = 'calling ' + str(fcn)
                [self.status, wa4] = self.call(fcn, self.params, functkw)
                if self.status < 0:
                    self.errmsg = 'WARNING: premature termination by "' + fcn + '"'
                    return
                fnorm1 = self.enorm(wa4)

                # Compute the scaled actual reduction
                catch_msg = 'computing convergence criteria'
                actred = -1.
                if (0.1 * fnorm1) < self.fnorm:
                    actred = - (fnorm1 / self.fnorm) ** 2 + 1.

                # Compute the scaled predicted reduction and the scaled directional
                # derivative
                for j in range(n):
                    wa3[j] = 0
                    wa3[0:j + 1] = wa3[0:j + 1] + fjac[0:j + 1, j] * wa1[
                        ipvt[j]]

                # Remember, alpha is the fraction of the full LM step actually
                # taken
                temp1 = self.enorm(alpha * wa3) / self.fnorm
                temp2 = (np.sqrt(alpha * par) * pnorm) / self.fnorm
                prered = temp1 * temp1 + (temp2 * temp2) / 0.5
                dirder = -(temp1 * temp1 + temp2 * temp2)

                # Compute the ratio of the actual to the predicted reduction.
                ratio = 0.
                if prered != 0:
                    ratio = actred / prered

                # Update the step bound
                if ratio <= 0.25:
                    if actred >= 0:
                        temp = .5
                    else:
                        temp = .5 * dirder / (dirder + .5 * actred)
                    if ((0.1 * fnorm1) >= self.fnorm) or (temp < 0.1):
                        temp = 0.1
                    delta = temp * np.min([delta, pnorm / 0.1])
                    par = par / temp
                else:
                    if (par == 0) or (ratio >= 0.75):
                        delta = pnorm / .5
                        par = .5 * par

                # Test for successful iteration
                if ratio >= 0.0001:
                    # Successful iteration.  Update x, fvec, and their norms
                    x = wa2
                    wa2 = diag * x
                    fvec = wa4
                    xnorm = self.enorm(wa2)
                    self.fnorm = fnorm1
                    self.niter = self.niter + 1

                # Tests for convergence
                if (np.abs(actred) <= ftol) and (prered <= ftol) \
                        and (0.5 * ratio <= 1):
                    self.status = 1
                if delta <= xtol * xnorm:
                    self.status = 2
                if (np.abs(actred) <= ftol) and (prered <= ftol) \
                        and (0.5 * ratio <= 1) and (self.status == 2):
                    self.status = 3
                if self.status != 0:
                    break

                # Tests for termination and stringent tolerances
                if self.niter >= maxiter:
                    self.status = 5
                if (np.abs(actred) <= machep) and (prered <= machep) \
                        and (0.5 * ratio <= 1):
                    self.status = 6
                if delta <= machep * xnorm:
                    self.status = 7
                if gnorm <= machep:
                    self.status = 8
                if self.status != 0:
                    break

                # End of inner loop. Repeat if iteration unsuccessful
                if ratio >= 0.0001:
                    break

                # Check for over/underflow
                if ~np.all(np.isfinite(wa1) & np.isfinite(wa2) & \
                           np.isfinite(x)) or ~np.isfinite(ratio):
                    errmsg = ('''ERROR: parameter or function value(s) have become
                        'infinite; check model function for over- 'and underflow''')
                    self.status = -16
                    break
                    # wh = where(finite(wa1) EQ 0 OR finite(wa2) EQ 0 OR finite(x) EQ 0, ct)
                    # if ct GT 0 OR finite(ratio) EQ 0 then begin

            if self.status != 0:
                break;
        # End of outer loop.

        catch_msg = 'in the termination phase'
        # Termination, either normal or user imposed.
        if len(self.params) == 0:
            return
        if nfree == 0:
            self.params = xall.copy()
        else:
            self.params[ifree] = x
        if (nprint > 0) and (self.status > 0):
            catch_msg = 'calling ' + str(fcn)
            [status, fvec] = self.call(fcn, self.params, functkw)
            catch_msg = 'in the termination phase'
            self.fnorm = self.enorm(fvec)

        if (self.fnorm is not None) and (fnorm1 is not None):
            self.fnorm = np.max([self.fnorm, fnorm1])
            self.fnorm = self.fnorm ** 2.

        self.covar = None
        self.perror = None
        # (very carefully) set the covariance matrix COVAR
        if (self.status > 0) and (nocovar == 0) and (n is not None) \
                and (fjac is not None) and (ipvt is not None):
            sz = fjac.shape
            if (n > 0) and (sz[0] >= n) and (sz[1] >= n) \
                    and (len(ipvt) >= n):

                catch_msg = 'computing the covariance matrix'
                cv = self.calc_covar(fjac[0:n, 0:n], ipvt[0:n])
                cv.shape = [n, n]
                nn = len(xall)

                # Fill in actual covariance matrix, accounting for fixed
                # parameters.
                self.covar = np.zeros([nn, nn], dtype=float)
                for i in range(n):
                    self.covar[ifree, ifree[i]] = cv[:, i]

                # Compute errors in parameters
                catch_msg = 'computing parameter errors'
                self.perror = np.zeros(nn, dtype=float)
                d = np.diagonal(self.covar).copy()
                wh = (np.nonzero(d >= 0))[0]
                if len(wh) > 0:
                    self.perror[wh] = np.sqrt(d[wh])
        return

    def __str__(self):
        return {'params': self.params,
                'niter': self.niter,
                'params': self.params,
                'covar': self.covar,
                'perror': self.perror,
                'status': self.status,
                'debug': self.debug,
                'errmsg': self.errmsg,
                'nfev': self.nfev,
                'damp': self.damp
                # ,'machar':self.machar
                }.__str__()

    # Default procedure to be called every iteration.  It simply prints
    # the parameter values.
    def defiter(self, fcn, x, iter, fnorm=None, functkw=None,
                quiet=0, iterstop=None, parinfo=None,
                format=None, pformat='%.10g', dof=1):
        """
        Print the current iteration information during optimization.

        Args:
            fcn (callable): The function to be minimized.
            x (array): Current parameter estimates.
            iter (int): The current iteration number.
            fnorm (float, optional): Current function norm (chi-squared).
                If not provided, it will be computed.
            functkw (dict, optional): Additional keyword arguments to be
                passed to `fcn`.
            quiet (int, optional): If set to a non-zero value, suppress output.
                Default is 0 (output shown).
            iterstop (optional): A hypothetical parameter for iteration stopping condition.
            parinfo (list of dicts, optional): Information for each parameter such as names.
            format (optional): Custom format for printing.
            pformat (str, optional): Format for printing parameters. Default is '%.10g'.
            dof (int, optional): Degrees of freedom. Default is 1.

        Returns:
            int: Always returns 0 after printing.
        """

        if self.debug:
            print()
            'Entering defiter...'
        if quiet:
            return
        if fnorm is None:
            [status, fvec] = self.call(fcn, x, functkw)
            fnorm = self.enorm(fvec) ** 2

        # Determine which parameters to print
        nprint = len(x)
        print()
        "Iter ", ('%6i' % iter), "   CHI-SQUARE = ", (
                '%.10g' % fnorm), " DOF = ", ('%i' % dof)
        for i in range(nprint):
            if (parinfo is not None) and ('parname' in parinfo[i]):
                p = '   ' + parinfo[i]['parname'] + ' = '
            else:
                p = '   P' + str(i) + ' = '
            if (parinfo is not None) and ('mpprint' in parinfo[i]):
                iprint = parinfo[i]['mpprint']
            else:
                iprint = 1
            if iprint:
                print()
                p + (pformat % x[i]) + '  '
        return 0

    #  DO_ITERSTOP:
    #  if keyword_set(iterstop) then begin
    #	  k = get_kbrd(0)
    #	  if k EQ string(byte(7)) then begin
    #		  message, 'WARNING: minimization not complete', /info
    #		  print, 'Do you want to terminate this procedure? (y/n)', $
    #			format='(A,$)'
    #		  k = ''
    #		  read, k
    #		  if strupcase(strmid(k,0,1)) EQ 'Y' then begin
    #			  message, 'WARNING: Procedure is terminating.', /info
    #			  mperr = -1
    #		  endif
    #	  endif
    #  endif

    # Procedure to parse the parameter values in PARINFO, which is a list of dictionaries
    def parinfo(self, parinfo=None, key='a', default=None, n=0):
        """
        Extract values from a parameter information structure.

        Args:
            parinfo (list of dicts, optional): List containing parameter information.
            key (str, optional): The key to look for within `parinfo`. Default is 'a'.
            default (optional): Default value to return if the key is not found.
            n (int, optional): The number of parameters to return. Defaults to 0.

        Returns:
            list: A list of values extracted from `parinfo`, or default values if not found.
        """

        if self.debug:
            print()
            'Entering parinfo...'
        if (n == 0) and (parinfo is not None):
            n = len(parinfo)
        if n == 0:
            values = default

            return values
        values = []
        for i in range(n):
            if (parinfo is not None) and (key in parinfo[i]):
                values.append(parinfo[i][key])
            else:
                values.append(default)

        # Convert to numeric arrays if possible
        test = default
        if type(default) == list:
            test = default[0]
        if isinstance(test, int):
            values = np.asarray(values, int)
        elif isinstance(test, float):
            values = np.asarray(values, float)
        return values

    # Call user function or procedure, with _EXTRA or not, with
    # derivatives or not.
    def call(self, fcn, x, functkw, fjac=None):
        """
        Call the user-defined function with the provided parameters.

        Args:
            fcn (callable): The user's function to be called.
            x (array): The current parameter estimates to be passed to `fcn`.
            functkw (dict): Additional keyword arguments to be passed to `fcn`.
            fjac (optional): Jacobian matrix or related parameter information.

        Returns:
            list: A list containing:
                - status (int): Status of the function call.
                - f (array): Function values calculated at `x`.
        """

        if self.debug:
            print()
            'Entering call...'
        if self.qanytied:
            x = self.tie(x, self.ptied)
        self.nfev = self.nfev + 1
        if fjac is None:
            [status, f] = fcn(x, fjac=fjac, **functkw)
            if self.damp > 0:
                # Apply the damping if requested.  This replaces the residuals
                # with their hyperbolic tangent.  Thus residuals larger than
                # DAMP are essentially clipped.
                f = np.tanh(f / self.damp)
            return [status, f]
        else:
            return fcn(x, fjac=fjac, **functkw)

    def enorm(self, vec):
        """
        Compute the Euclidean norm of a vector.

        Args:
            vec (array): The input vector.

        Returns:
            float: The Euclidean norm of the input vector.
        """
        ans = self.blas_enorm(vec)
        return ans

    def fdjac2(self, fcn, x, fvec, step=None, ulimited=None, ulimit=None,
               dside=None, epsfcn=None, autoderivative=1,
               functkw=None, xall=None, ifree=None, dstep=None):
        """
        Compute the Jacobian matrix using finite differences.

        Args:
            fcn (callable): The function to compute the Jacobian for.
            x (array): The current parameter estimates.
            fvec (array): Function values at `x`.
            step (array, optional): Steps for finite differencing. Default is None.
            ulimited (array, optional): Flags indicating whether parameters are unlimited.
            ulimit (array, optional): Upper limits for parameters.
            dside (array, optional): Side of finite difference to use.
            epsfcn (float, optional): Parameter for finite differencing.
            autoderivative (int, optional): Set to 0 to use analytical derivatives instead.
            functkw (dict, optional): Additional keyword arguments for `fcn`.
            xall (array, optional): Array of all parameters.
            ifree (array, optional): Indices of free parameters.
            dstep (array, optional): Absolute or relative step size.

        Returns:
            array: Jacobian matrix evaluated at `x`.
        """
        if self.debug:
            print()
            'Entering fdjac2...'
        machep = self.machar.machep
        if epsfcn is None:
            epsfcn = machep
        if xall is None:
            xall = x
        if ifree is None:
            ifree = np.arange(len(xall))
        if step is None:
            step = x * 0.
        nall = len(xall)

        eps = np.sqrt(np.max([epsfcn, machep]))
        m = len(fvec)
        n = len(x)

        # Compute analytical derivative if requested
        if autoderivative == 0:
            mperr = 0
            fjac = np.zeros(nall, dtype=float)
            fjac[ifree] = 1.0  # Specify which parameters need derivatives
            [status, fp] = self.call(fcn, xall, functkw, fjac=fjac)

            if len(fjac) != m * nall:
                print()
                'ERROR: Derivative matrix was not computed properly.'
                return None

            # This definition is consistent with CURVEFIT
            # Sign error found (thanks Jesus Fernandez <fernande@irm.chu-caen.fr>)
            fjac.shape = [m, nall]
            fjac = -fjac

            # Select only the free parameters
            if len(ifree) < nall:
                fjac = fjac[:, ifree]
                fjac.shape = [m, n]
                return fjac

        fjac = np.zeros([m, n], dtype=float)

        h = eps * np.abs(x)

        # if STEP is given, use that
        # STEP includes the fixed parameters
        if step is not None:
            stepi = step[ifree]
            wh = (np.nonzero(stepi > 0))[0]
            if len(wh) > 0:
                h[wh] = stepi[wh]

        # if relative step is given, use that
        # DSTEP includes the fixed parameters
        if len(dstep) > 0:
            dstepi = dstep[ifree]
            wh = (np.nonzero(dstepi > 0))[0]
            if len(wh) > 0:
                h[wh] = np.abs(dstepi[wh] * x[wh])

        # In case any of the step values are zero
        h[h == 0] = eps

        # Reverse the sign of the step if we are up against the parameter
        # limit, or if the user requested it.
        # DSIDE includes the fixed parameters (ULIMITED/ULIMIT have only
        # varying ones)
        mask = dside[ifree] == -1
        if len(ulimited) > 0 and len(ulimit) > 0:
            mask = (mask | ((ulimited != 0) & (x > ulimit - h)))
            wh = (np.nonzero(mask))[0]
            if len(wh) > 0:
                h[wh] = - h[wh]
        # Loop through parameters, computing the derivative for each
        for j in range(n):
            xp = xall.copy()
            xp[ifree[j]] = xp[ifree[j]] + h[j]
            [status, fp] = self.call(fcn, xp, functkw)
            if status < 0:
                return None

            if np.abs(dside[ifree[j]]) <= 1:
                # COMPUTE THE ONE-SIDED DERIVATIVE
                # Note optimization fjac(0:*,j)
                fjac[0:, j] = (fp - fvec) / h[j]

            else:
                # COMPUTE THE TWO-SIDED DERIVATIVE
                xp[ifree[j]] = xall[ifree[j]] - h[j]

                mperr = 0
                [status, fm] = self.call(fcn, xp, functkw)
                if status < 0:
                    return None

                # Note optimization fjac(0:*,j)
                fjac[0:, j] = (fp - fm) / (2 * h[j])
        return fjac

    def qrfac(self, a, pivot=0):
        """
        Perform QR factorization of a matrix.

        Args:
            a (array): Input matrix to be factorized.
            pivot (int, optional): Whether to perform pivoting (0 for no, 1 for yes).

        Returns:
            list: A list containing:
                - a (array): The upper triangular matrix after factorization.
                - ipvt (array): The pivot indices.
                - rdiag (array): The diagonal elements of the R matrix.
                - acnorm (array): The norms of the columns of A.
        """

        if self.debug: print()
        'Entering qrfac...'
        machep = self.machar.machep
        sz = a.shape
        m = sz[0]
        n = sz[1]

        # Compute the initial column norms and initialize arrays
        acnorm = np.zeros(n, dtype=float)
        for j in range(n):
            # CALCULATE THE NORM OF EACH COLUMN IN JACOBIAN
            # STORE IN ACNORM 1Xn array
            acnorm[j] = self.enorm(a[:, j])
        rdiag = acnorm.copy()
        wa = rdiag.copy()
        ipvt = np.arange(n)

        # Reduce a to r with householder transformations
        minmn = np.min([m, n])
        for j in range(minmn):
            if pivot != 0:
                # Bring the column of largest norm into the pivot position
                rmax = np.max(rdiag[j:])
                kmax = (np.nonzero(rdiag[j:] == rmax))[0]
                ct = len(kmax)
                kmax = kmax + j
                if ct > 0:
                    kmax = kmax[0]

                    # Exchange rows via the pivot only.  Avoid actually exchanging
                    # the rows, in case there is lots of memory transfer.  The
                    # exchange occurs later, within the body of MPFIT, after the
                    # extraneous columns of the matrix have been shed.
                    if kmax != j:
                        temp = ipvt[j];
                        ipvt[j] = ipvt[kmax];
                        ipvt[kmax] = temp
                        rdiag[kmax] = rdiag[j]
                        wa[kmax] = wa[j]

            # Compute the householder transformation to reduce the jth
            # column of A to a multiple of the jth unit vector
            lj = ipvt[j]
            ajj = a[j:, lj]
            ajnorm = self.enorm(ajj)
            if ajnorm == 0:
                break
            if a[j, lj] < 0:
                ajnorm = -ajnorm

            ajj = ajj / ajnorm
            ajj[0] = ajj[0] + 1
            # *** Note optimization a(j:*,j)
            a[j:, lj] = ajj

            # Apply the transformation to the remaining columns
            # and update the norms

            # NOTE to SELF: tried to optimize this by removing the loop,
            # but it actually got slower.  Reverted to "for" loop to keep
            # it simple.
            if j + 1 < n:
                for k in range(j + 1, n):
                    lk = ipvt[k]
                    ajk = a[j:, lk]
                    # *** Note optimization a(j:*,lk)
                    # (corrected 20 Jul 2000)
                    if a[j, lj] != 0:
                        # CALCULATING GRAM-SCHMIDT
                        a[j:, lk] = ajk - ajj * sum(ajk * ajj) / a[j, lj]
                        if (pivot != 0) and (rdiag[k] != 0):
                            temp = a[j, lk] / rdiag[k]
                            rdiag[k] = rdiag[k] * np.sqrt(
                                np.max([(1. - temp ** 2), 0.]))
                            temp = rdiag[k] / wa[k]
                            if (0.05 * temp * temp) <= machep:
                                rdiag[k] = self.enorm(a[j + 1:, lk])
                                wa[k] = rdiag[k]
            rdiag[j] = -ajnorm
        return [a, ipvt, rdiag, acnorm]

    def qrsolv(self, r, ipvt, diag, qtb, sdiag):
        """
        Solve a linear system using QR factorization.

        Args:
            r (array): Upper triangular matrix from QR factorization.
            ipvt (array): Pivot indices.
            diag (array): Diagonal scaling factors.
            qtb (array): The product of (Q^T) and the right-hand side vector b.
            sdiag (array): Storage for diagonal elements.

        Returns:
            tuple: A tuple containing:
                - r (array): Modified upper triangular matrix.
                - x (array): Solution vector.
                - sdiag (array): Updated diagonal elements.
        """

        if self.debug:
            print()
            'Entering qrsolv...'
        sz = r.shape
        m = sz[0]
        n = sz[1]

        # copy r and (q transpose)*b to preserve input and initialize s.
        # in particular, save the diagonal elements of r in x.

        for j in range(n):
            r[j:n, j] = r[j, j:n]
        x = np.diagonal(r).copy()
        wa = qtb.copy()

        # Eliminate the diagonal matrix d using a givens rotation
        for j in range(n):
            l = ipvt[j]
            if diag[l] == 0:
                break
            sdiag[j:] = 0
            sdiag[j] = diag[l]

            # The transformations to eliminate the row of d modify only a
            # single element of (q transpose)*b beyond the first n, which
            # is initially zero.

            qtbpj = 0.
            for k in range(j, n):
                if sdiag[k] == 0:
                    break
                if np.abs(r[k, k]) < np.abs(sdiag[k]):
                    cotan = r[k, k] / sdiag[k]
                    sine = 0.5 / np.sqrt(.25 + .25 * cotan * cotan)
                    cosine = sine * cotan
                else:
                    tang = sdiag[k] / r[k, k]
                    cosine = 0.5 / np.sqrt(.25 + .25 * tang * tang)
                    sine = cosine * tang

                # Compute the modified diagonal element of r and the
                # modified element of ((q transpose)*b,0).
                r[k, k] = cosine * r[k, k] + sine * sdiag[k]
                temp = cosine * wa[k] + sine * qtbpj
                qtbpj = -sine * wa[k] + cosine * qtbpj
                wa[k] = temp

                # Accumulate the transformation in the row of s
                if n > k + 1:
                    temp = cosine * r[k + 1:n, k] + sine * sdiag[k + 1:n]
                    sdiag[k + 1:n] = -sine * r[k + 1:n, k] + cosine * sdiag[
                                                                      k + 1:n]
                    r[k + 1:n, k] = temp
            sdiag[j] = r[j, j]
            r[j, j] = x[j]

        # Solve the triangular system for z.  If the system is singular
        # then obtain a least squares solution
        nsing = n
        wh = (np.nonzero(sdiag == 0))[0]
        if len(wh) > 0:
            nsing = wh[0]
            wa[nsing:] = 0

        if nsing >= 1:
            wa[nsing - 1] = wa[nsing - 1] / sdiag[nsing - 1]  # Degenerate case
            # *** Reverse loop ***
            for j in range(nsing - 2, -1, -1):
                sum0 = sum(r[j + 1:nsing, j] * wa[j + 1:nsing])
                wa[j] = (wa[j] - sum0) / sdiag[j]

        # Permute the components of z back to components of x
        x[ipvt] = wa
        return (r, x, sdiag)

    def lmpar(self, r, ipvt, diag, qtb, delta, x, sdiag, par=None):
        """
        Perform damped least squares minimization.

        Args:
            r (array): Upper triangular matrix from QR factorization.
            ipvt (array): Pivot indices.
            diag (array): Diagonal scaling factors.
            qtb (array): The product of (Q^T) and the right-hand side vector b.
            delta (float): Convergence threshold.
            x (array): Initial guess for the solution.
            sdiag (array): Storage for diagonal elements from the QR factorization.
            par (float, optional): Initial value for the parameter.

        Returns:
            list: A list containing:
                - r (array): Modified upper triangular matrix.
                - par (float): Updated parameter value.
                - x (array): Solution vector.
                - sdiag (array): Updated diagonal elements.
        """

        if self.debug:
            print()
            'Entering lmpar...'
        dwarf = self.machar.minnum
        machep = self.machar.machep
        sz = r.shape
        m = sz[0]
        n = sz[1]

        # Compute and store in x the gauss-newton direction.  If the
        # jacobian is rank-deficient, obtain a least-squares solution
        nsing = n
        wa1 = qtb.copy()
        rthresh = np.max(np.abs(np.diagonal(r).copy())) * machep
        wh = (np.nonzero(np.abs(np.diagonal(r).copy()) < rthresh))[0]
        if len(wh) > 0:
            nsing = wh[0]
            wa1[wh[0]:] = 0
        if nsing >= 1:
            # *** Reverse loop ***
            for j in range(nsing - 1, -1, -1):
                wa1[j] = wa1[j] / r[j, j]
                if j - 1 >= 0:
                    wa1[0:j] = wa1[0:j] - r[0:j, j] * wa1[j]

        # Note: ipvt here is a permutation array
        x[ipvt] = wa1

        # Initialize the iteration counter.  Evaluate the function at the
        # origin, and test for acceptance of the gauss-newton direction
        iter = 0
        wa2 = diag * x
        dxnorm = self.enorm(wa2)
        fp = dxnorm - delta
        if fp <= 0.1 * delta:
            return [r, 0., x, sdiag]

        # If the jacobian is not rank deficient, the newton step provides a
        # lower bound, parl, for the zero of the function.  Otherwise set
        # this bound to zero.

        parl = 0.
        if nsing >= n:
            wa1 = diag[ipvt] * wa2[ipvt] / dxnorm
            wa1[0] = wa1[0] / r[0, 0]  # Degenerate case
            for j in range(1, n):  # Note "1" here, not zero
                sum0 = sum(r[0:j, j] * wa1[0:j])
                wa1[j] = (wa1[j] - sum0) / r[j, j]

            temp = self.enorm(wa1)
            parl = ((fp / delta) / temp) / temp

        # Calculate an upper bound, paru, for the zero of the function
        for j in range(n):
            sum0 = sum(r[0:j + 1, j] * qtb[0:j + 1])
            wa1[j] = sum0 / diag[ipvt[j]]
        gnorm = self.enorm(wa1)
        paru = gnorm / delta
        if paru == 0:
            paru = dwarf / np.min([delta, 0.1])

        # If the input par lies outside of the interval (parl,paru), set
        # par to the closer endpoint

        par = np.max([par, parl])
        par = np.min([par, paru])
        if par == 0:
            par = gnorm / dxnorm

        # Beginning of an interation
        while (1):
            iter = iter + 1

            # Evaluate the function at the current value of par
            if par == 0:
                par = np.max([dwarf, paru * 0.001])
            temp = np.sqrt(par)
            wa1 = temp * diag
            [r, x, sdiag] = self.qrsolv(r, ipvt, wa1, qtb, sdiag)
            wa2 = diag * x
            dxnorm = self.enorm(wa2)
            temp = fp
            fp = dxnorm - delta

            if (np.abs(fp) <= 0.1 * delta) or \
                    ((parl == 0) and (fp <= temp) and (temp < 0)) or \
                    (iter == 10):
                break;

            # Compute the newton correction
            wa1 = diag[ipvt] * wa2[ipvt] / dxnorm

            for j in range(n - 1):
                wa1[j] = wa1[j] / sdiag[j]
                wa1[j + 1:n] = wa1[j + 1:n] - r[j + 1:n, j] * wa1[j]
            wa1[n - 1] = wa1[n - 1] / sdiag[n - 1]  # Degenerate case

            temp = self.enorm(wa1)
            parc = ((fp / delta) / temp) / temp

            # Depending on the sign of the function, update parl or paru
            if fp > 0:
                parl = np.max([parl, par])
            if fp < 0:
                paru = np.min([paru, par])

            # Compute an improved estimate for par
            par = np.max([parl, par + parc])

            # End of an iteration

        # Termination
        return [r, par, x, sdiag]

    # Procedure to tie one parameter to another.
    def tie(self, p, ptied=None):
        """
        Apply constraints by tying parameters to each other.

        Args:
            p (array): Current parameters.
            ptied (list of str, optional): List of parameter names to tie.

        Returns:
            array: The modified parameters.
        """

        if self.debug:
            print()
            'Entering tie...'
        if ptied is None:
            return
        for i in range(len(ptied)):
            if ptied[i] == '':
                continue
            cmd = 'p[' + str(i) + '] = ' + ptied[i]
            exec(cmd)
        return p

    def calc_covar(self, rr, ipvt=None, tol=1.e-14):
        """
        Calculate the covariance matrix based on the QR factorization.

        Args:
            rr (array): The upper triangular matrix from QR factorization.
            ipvt (array, optional): Pivot indices.
            tol (float, optional): Tolerance for determining if a parameter is significant.

        Returns:
            array: The covariance matrix or -1 in case of an error.
        """

        if self.debug:
            print()
            'Entering calc_covar...'
        if rr.ndim != 2:
            print()
            'ERROR: r must be a two-dimensional matrix'
            return -1
        s = rr.shape
        n = s[0]
        if s[0] != s[1]:
            print()
            'ERROR: r must be a square matrix'
            return -1

        if ipvt is None:
            ipvt = np.arange(n)
        r = rr.copy()
        r.shape = [n, n]

        # For the inverse of r in the full upper triangle of r
        l = -1
        tolr = tol * np.abs(r[0, 0])
        for k in range(n):
            if np.abs(r[k, k]) <= tolr:
                break
            r[k, k] = 1. / r[k, k]
            for j in range(k):
                temp = r[k, k] * r[j, k]
                r[j, k] = 0.
                r[0:j + 1, k] = r[0:j + 1, k] - temp * r[0:j + 1, j]
            l = k

        # Form the full upper triangle of the inverse of (r transpose)*r
        # in the full upper triangle of r
        if l >= 0:
            for k in range(l + 1):
                for j in range(k):
                    temp = r[j, k]
                    r[0:j + 1, j] = r[0:j + 1, j] + temp * r[0:j + 1, k]
                temp = r[k, k]
                r[0:k + 1, k] = temp * r[0:k + 1, k]

        # For the full lower triangle of the covariance matrix
        # in the strict lower triangle or and in wa
        wa = np.repeat([r[0, 0]], n)
        for j in range(n):
            jj = ipvt[j]
            sing = j > l
            for i in range(j + 1):
                if sing:
                    r[i, j] = 0.
                ii = ipvt[i]
                if ii > jj:
                    r[ii, jj] = r[i, j]
                if ii < jj:
                    r[jj, ii] = r[i, j]
            wa[jj] = r[j, j]

        # Symmetrize the covariance matrix in r
        for j in range(n):
            r[0:j + 1, j] = r[j, 0:j + 1]
            r[j, j] = wa[j]

        return r


class machar:
    """
    A class to hold machine learning variables related to floating point numbers.

    This class gathers important constants and properties of floating point
    representation, specifically for 32-bit and 64-bit float types,
    using NumPy's `finfo` to provide values that characterize the limits
    of floating point precision.

    Attributes:
        machep (float): The machine epsilon, the smallest value such that
            `1.0 + machep != 1.0` for the specified floating point type.
        maxnum (float): The maximum representable positive floating point
            number for the specified type.
        minnum (float): The minimum representable positive floating point
            number (i.e., the smallest positive number greater than zero).
        maxlog (float): The natural logarithm of the maximum representable
            positive floating point number.
        minlog (float): The natural logarithm of the minimum representable
            positive floating point number.
        rdwarf (float): A parameter used in scaling small numbers, defined as
            `sqrt(minnum * 1.5) * 10`.
        rgiant (float): A parameter used in scaling large numbers, defined as
            `sqrt(maxnum) * 0.1`.

    Args:
        double (int, optional): An indicator for the floating point type to
            use. If set to 0, 32-bit floats (`float32`) are utilized;
            if set to 1 (or any other value), 64-bit floats (`float64`) are used.
            Default is 1 (using 64-bit floats).

    Example:
        >>> m = machar(double=1)
        >>> print(m.machep)
        2.220446049250313e-16
    """

    def __init__(self, double=1):
        if double == 0:
            info = np.finfo(np.float32)
        else:
            info = np.finfo(np.float64)

        self.machep = info.eps
        self.maxnum = info.max
        self.minnum = info.tiny

        self.maxlog = np.log(self.maxnum)
        self.minlog = np.log(self.minnum)
        self.rdwarf = np.sqrt(self.minnum * 1.5) * 10
        self.rgiant = np.sqrt(self.maxnum) * 0.1
