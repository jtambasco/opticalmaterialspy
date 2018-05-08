import numpy as np
from scipy import constants as spc
from scipy import interpolate as spi
import abc


class _Material(metaclass=abc.ABCMeta):
    def __init__(self, wlMinNm=300., wlMaxNm=5000.):
        self._wlMin = wlMinNm
        self._wlMax = wlMaxNm

        self._wls = np.linspace(self._wlMin, self._wlMax, 100, False)
        self._wlStep = (self._wlMax-self._wlMin)/100.

        self.nDer1Func = None
        self.nDer2Func = None

    def _convertWavelengthUnitsNm(self, wavelength):
        if wavelength is not False and wavelength is not None:
            if np.all(wavelength <= self._wlMax * 1.e-9):
                wlFactor = 1.e9
            elif np.all(wavelength <= self._wlMax * 1.e-3):
                wlFactor = 1.e3
            else:
                wlFactor = 1.
            wl = wavelength * wlFactor
            assert(np.all(wl >= self._wlMin) and np.all(wl <= self._wlMax))
        else:
            wl = wavelength
        return wl

    def convertWavelengthUnitsNm(func):
        def _convertWavelengthUnitsNmWrapper(self, wavelength):
            wl = self._convertWavelengthUnitsNm(wavelength)
            return func(self, wl)
        return _convertWavelengthUnitsNmWrapper

    @abc.abstractproperty
    def _eps(self, wavelength=None):
        return NotImplementedError

    @convertWavelengthUnitsNm
    def eps(self, wavelength=None):
        return self._eps(wavelength)

    def _nDer1(self, wavelength):
        if not self.nDer1Func:
            nDer1 = np.gradient(self.n(self._wls), self._wlStep) # [1/nm]
            self.nDer1Func = spi.interp1d(self._wls, nDer1)
        return self.nDer1Func(wavelength) * 1.e9  # [1/nm] -> [1/m]

    def _nDer2(self, wavelength):
        if not self.nDer2Func:
            nDer2 = np.gradient(self._nDer1(self._wls), self._wlStep)  # [1/(nm)^2]
            self.nDer2Func = spi.interp1d(self._wls, nDer2)
        return self.nDer2Func(wavelength) * 1.e9  # [1/(nm)^2] -> [1/m]

    def n(self, wavelength=None):
        return self.eps(wavelength)**0.5

    # dn/dlambda
    def nDer1(self, wavelength):
        return self._nDer1(wavelength)

    # d^2n/dlambda^2
    def nDer2(self, wavelength):
        return self._nDer2(wavelength)

    # n_g: group velocity index
    @convertWavelengthUnitsNm
    def ng(self, wavelength):
        return self.n(wavelength) - (wavelength*1.e-9)*self.nDer1(wavelength)

    # Group velocity (v_g)
    def vg(self, wavelength):
        return spc.c / self.ng(wavelength)

    # Group velocity dispersion (GVD)
    # beta_2 = lambda^3/(2 pi c^2) * d^2n/dlambda^2
    @convertWavelengthUnitsNm
    def gvd(self, wavelength):
        g = (wavelength*1.e-9)**3./(2.*spc.pi*spc.c**2.) * self.nDer2(wavelength)
        return g

    # beta(w) = beta0 + beta1(w-w0) + 1/2 beta2(w-w0)^2 + ...
    @convertWavelengthUnitsNm
    def beta0(self, wavelength):
        return 2.*spc.pi*self.n(wavelength)/(wavelength*1.e-9)

    def beta1(self, wavelength):
        return 1./self.vg(wavelength)

    def beta2(self, wavelength):
        return self.gvd(wavelength)

    @staticmethod
    def cauchy_equation(wavelength, coeffients):
        n = 0.
        for i, c in enumerate(coeffients):
            exponent  = 2*i
            n += c / wavelength**exponent
        return n

class _MaterialAni(metaclass=abc.ABCMeta):
    def __init__(self):
        pass

    @abc.abstractproperty
    def xx(self):
        pass

    @abc.abstractproperty
    def yy(self):
        pass

    @abc.abstractproperty
    def zz(self):
        pass

    @property
    def xy(self):
        return None

    @property
    def xz(self):
        return None

    @property
    def yx(self):
        return None

    @property
    def yz(self):
        return None

    @property
    def zx(self):
        return None

    @property
    def zy(self):
        return None
