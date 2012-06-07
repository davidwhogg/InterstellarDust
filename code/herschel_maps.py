"""
This code is part of the Interstellar Dust project.

Copyright 2012 David W. Hogg

This code builds a simultaneous emitting dust model for a set of
aligned Herschel maps.

It has various things related to Herschel data hard-coded.
"""

import numpy as np
from scipy.signal import fftconvolve as convolve
import scipy.optimize as op
import pyfits as pyf

# all the following from physics.nist.gov
hh = 6.62606957e-34 # J s
cc = 299792458. # m s^{-1}
kk = 1.3806488e-23 # J K^{-1}

# MAGIC number: wavelength scale for emissivity model
lam0 = 1.e-4 # m

def black_body(lam, lnT):
    """
    Compute the black-body formula, for a given lnT.
    """
    return (2. * hh * cc ** 2 * lam ** -5 /
            (np.exp(hh * cc / (lam * kk * np.exp(lnT))) - 1.))

def synthesize_map(lam, lnRho, lnT, beta):
    """
    Produce a synthetic Herschel map.

    ### issues:
    - assumes maps are single-wavelength (narrow-band)
    - assumes there is only one temperature component per pixel
    - has an idiotic emissivity model
    - kernel convolution is slow-ass
    - units for lambda, rho, T all implicit

    ### inputs:
    - `lam`: central wavelength for the map
    - `lnRho`: two-d map of (log) dust density
    - `lnT`: scalar or else two-d map of (log) dust temperature
    - `beta`: scalar or else two-d map of emissivity parameter

    ### outputs:
    - `map`: a two-d map of (log) intensity, same size as lnRho
    """
    return lnRho + np.log(black_body(lam, lnT)) + np.log(lam / lam0) * beta

def chi(data, pars, prior_info):
    """
    Produce a single-image residual map appropriate for least-square
    fitting.
    """
    chi = np.array([])
    data.set_calibration(pars.get_calibration_pars())
    for image, invvar, sky, lam, kernel in data:
        thischi = (np.sqrt(invvar) *
                   (image - sky - convolve(pars.get_map(lam), kernel, mode="valid")))
        chi = np.append(chi, thischi)
    return chi

def resampling_xy(shape, WCS, refWCS):
    """
    Get pixel positions in refWCS corresponding to the pixels in WCS.
    """
    x = np.arange(shape[0])[:, None]
    y = np.arange(shape[1])[None, :]
    refx, refy = refWCS.position2pixel(WCS.pixel2position(x, y))
    return (refx, refy)

class Data():

    def __init__(self, dir):
        """
        Brittle function to read Groves data sample and make the
        things we need for fitting.

        # issues:
        - Need to read the WCS too.
        - Need, for each image, to construct the rectangular
          nearest-neighbor interpolation indices to image 0.
        - Can I just blow up the SPIRE images with interpolation?
        """
        self.images = []
        self.lams = []
        self.WCSs = []
        dataList = [
            ('m31_brick15_PACS100.fits', 100.), # mum
            ('m31_brick15_PACS160.fits', 160.),
            ('m31_brick15_SPIRE250.fits', 250.),
            ('m31_brick15_SPIRE350.fits', 350.),
            ('m31_brick15_SPIRE500.fits', 500.),
            ]
        for n, (fn, lam) in enumerate(dataList):
            f = pyf.open(dir + fn)
            thisimage = np.array(f[0].data, dtype=float)
            thisWCS = get_wcs(f) # replace this with something NOT WRONG
            self.WCS = self.WCS + [thisWCS]
            f.close()
            self.images = self.images + [thisimage, ]
            self.lams = self.lams + [lam * 1.e-6, ] # m
            self.WCS = self.WCS + [thisWCS]
            self.xylist = resampling_xy(thisimage.shape, thisWCS, self.images[0].shape, self.WCS[0])
            print thisimage.shape, len(self.images), len(self.lams)
            print np.min(thisimage), np.median(thisimage), np.max(thisimage)
        return None

    def __getitem__(self, n):
        """
        Return what the optimizer needs from the data.
        """
        return self.images[n], self.invvars[n], self.skies[n], self.lams[n], self.kernels[n]

    def set_calibration(self, skyVariances, signalVariances, kernelVariances):
        """
        Because the calibration parameters relating to noise and PSF
        are fitting parameters, we need this function to update the
        current beliefs about calibration.

        This code has hard-coded the fact that images[0] is the
        highest resolution band.
        """
        for n in len(self.lams):
            self.invvar[n] = 1. / (skyVariances[n] + signalVariances[n] * self.images[n])
            hw = 9
            if k == 0:
                thiskernel = np.zeros((hw + hw + 1, hw + hw + 1))
                thiskernel[hw, hw] = 1.
            else:
                thiskernel = np.exp(-0.5 * (np.arange(-hw, hw+1)[:, None] ** 2 +
                                            np.arange(-hw, hw+1)[None, :] ** 2) /
                                     self.kernelVariances[n])
                thiskernel /= np.sum(thiskernel)
            self.kernel[n] = thiskernel
        pass

class Pars():

    def __init__(self, lnRho, lnT, beta, skyLevels, skyVariances, signalVariances, kernelVariances):
        self.lnRho = lnRho
        self.lnT = lnT
        self.beta = beta
        self.skyLevels = skyLevels
        self.skyVariances = skyVariances
        self.signalVariances = signalVariances
        self.kernelVariances = kernelVariances
        return None

    def get_map_pars(self):
        return self.lnRho, self.lnT, self.beta

    def get_map(self, lam):
        return synthesize_map(lam, *self.get_map_pars())

    def get_calibration_pars(self):
        return self.skyVariances, self.signalVariances, self.kernelVariances

    def get_lnRho_pars(self):
        return self.lnRho.flatten

    def set_lnRho_pars(self, pars):
        self.lnRho = pars.reshape(self.lnRho.shape)
        return None

    def get_kernel_pars(self):
        return self.kernelVariances

    def set_kernel_pars(self, pars):
        self.kernelVariances = pars
        return None

if __name__ == '__main__':
    data = Data('../data/')

