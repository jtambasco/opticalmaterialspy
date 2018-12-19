import numpy as np
from scipy import constants as spc
from scipy import interpolate as spi
import abc
from functools import wraps


class _Material(metaclass=abc.ABCMeta):
    '''
    Abstract material class that can calculate many
    material properties based on the permittivity.

    Args:
        wlMinNm (float): The minimum wavelength [nm].
        wlMaxNm (float): The maximum wavelength [nm].
    '''
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
            assert(np.all(wl >= self._wlMin) and np.all(wl <= self._wlMax)), \
                'Wavelength not in range %.1f nm to %.1f nm.' % (self._wlMin, self._wlMax)
        else:
            wl = wavelength
        return wl

    def convertWavelengthUnitsNm(func):
        @wraps(func)
        def _convertWavelengthUnitsNmWrapper(self, wavelength):
            wl = self._convertWavelengthUnitsNm(wavelength)
            return func(self, wl)
        return _convertWavelengthUnitsNmWrapper

    @abc.abstractproperty
    def _eps(self, wavelength=None):
        '''
        The permittivty of the desired material.

        Each new material needs this method implemented.

        Args:
            wavelength (float, list, None): The wavelength the permittivty
                will be evaluated at.

        Returns:
            float, list: The permittivty at the target wavelength.
        '''
        return NotImplementedError

    @convertWavelengthUnitsNm
    def eps(self, wavelength=None):
        '''
        The permittivty of the desired material.

        Args:
            wavelength (float, list, None): The wavelength the permittivty
                will be evaluated at.

        Returns:
            float, list: The permittivty at the target wavelength.
        '''
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
        '''
        The refractive index of the desired material.

        Args:
            wavelength (float, list, None): The wavelength the refractive index
                will be evaluated at.

        Returns:
            float, list: The refractive index at the target wavelength.
        '''
        return self.eps(wavelength)**0.5

    # dn/dlambda
    def nDer1(self, wavelength):
        '''
        The first derivative of the refractive index with respect to
        wavelength.

        Args:
            wavelength (float, list, None): The wavelength(s) the derivative
                will be evaluated at.

        Returns:
            float, list: The refractive index at the target wavelength(s).
        '''
        return self._nDer1(wavelength)

    # d^2n/dlambda^2
    def nDer2(self, wavelength):
        '''
        The second derivative of the refractive index with respect to
        wavelength.

        Args:
            wavelength (float, list, None): The wavelength(s) the derivative
                will be evaluated at.

        Returns:
            float, list: The refractive index at the target wavelength(s).
        '''
        return self._nDer2(wavelength)

    # n_g: group velocity index
    @convertWavelengthUnitsNm
    def ng(self, wavelength):
        '''
        The group index with respect to wavelength.

        Args:
            wavelength (float, list, None): The wavelength(s) the group
                index will be evaluated at.

        Returns:
            float, list: The group index at the target wavelength(s).
        '''
        return self.n(wavelength) - (wavelength*1.e-9)*self.nDer1(wavelength)

    # Group velocity (v_g)
    def vg(self, wavelength):
        '''
        The group velocities with respect to wavelength.

        Args:
            wavelength (float, list, None): The wavelength(s) the group
                velocities will be evaluated at.

        Returns:
            float, list: The group velocities at the target wavelength(s).
        '''
        return spc.c / self.ng(wavelength)

    # Group velocity dispersion (GVD)
    # beta_2 = lambda^3/(2 pi c^2) * d^2n/dlambda^2
    @convertWavelengthUnitsNm
    def gvd(self, wavelength):
        '''
        The group velocity dispersion (GVD) with respect to wavelength.

        Args:
            wavelength (float, list, None): The wavelength(s) the GVD will
                be evaluated at.

        Returns:
            float, list: The GVD at the target wavelength(s).
        '''
        g = (wavelength*1.e-9)**3./(2.*spc.pi*spc.c**2.) * self.nDer2(wavelength)
        return g

    # beta(w) = beta0 + beta1(w-w0) + 1/2 beta2(w-w0)^2 + ...
    @convertWavelengthUnitsNm
    def beta0(self, wavelength):
        '''
        The propagation constant with respect to wavelength.

        Args:
            wavelength (float, list, None): The wavelength(s) the
                propagation constant will be evaluated at.

        Returns:
            float, list: The propagation constant at the target wavelength(s).
        '''
        return 2.*spc.pi*self.n(wavelength)/(wavelength*1.e-9)

    def beta1(self, wavelength):
        '''
        The derivative of the propagation constant with respect to wavelength.

        Args:
            wavelength (float, list, None): The wavelength(s) the
                propagation constant will be evaluated at.

        Returns:
            float, list: The propagation constant at the target wavelength(s).
        '''
        return 1./self.vg(wavelength)

    def beta2(self, wavelength):
        '''
        The second derivative of the propagation constant with respect to wavelength.

        Args:
            wavelength (float, list, None): The wavelength(s) the
                propagation constant will be evaluated at.

        Returns:
            float, list: The propagation constant at the target wavelength(s).
        '''
        return self.gvd(wavelength)

    def z0(self, wavelength):
        '''
        The wave impedance assuming the material is dielectric (not
        lossy or magnetic).

        Args:
            wavelength (float, list, None): The wavelength(s) the
                propagation constant will be evaluated at.

        Returns:
            float, list: The impedance of the material.

        '''
        return 120*np.pi / self.n(wavelength)

    @staticmethod
    def _cauchy_equation(wavelength, coefficients):
        '''
        Helpful function to evaluate Cauchy equations.

        Args:
            wavelength (float, list, None): The wavelength(s) the
                Cauchy equation will be evaluated at.
            coefficients (list): A list of the coefficients of
                the Cauchy equation.

        Returns:
            float, list: The refractive index at the target wavelength(s).
        '''
        n = 0.
        for i, c in enumerate(coefficients):
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

    def n3(self, wl):
        return [self.xx.n(wl), self.yy.n(wl), self.zz.n(wl)]

    def n_xyz(self, wl):
        return self.n3(wl)

    def n5(self, wl):
        return [self.xx.n(wl), self.xy.n(wl), self.yx.n(wl),
                self.yy.n(wl), self.zz.n(wl)]

