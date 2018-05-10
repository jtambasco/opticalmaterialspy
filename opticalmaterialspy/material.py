# -*- coding: utf-8 -*-

import numpy as np
from scipy import constants as spc
from scipy import interpolate as spi
from ._material_base import _Material


class File(_Material):
    '''  `filename` is the name of a file with column 1 as the wavelength [nm],
         column 2 as the effective index of mode 1, column 3 as the effective
         index of mode 2 etc. '''
    def __init__(self, filename, mode=0, delimiter=','):
        txt = np.loadtxt(filename, delimiter=delimiter)
        wls = txt[:,0] * 1.e-9
        nEffs = txt[:,mode+1]

        wlMin = wls[0]
        wlMax = wls[-1]
        nEffFunc = spi.interp1d(wls, nEffs)
        self._eps = lambda wavelength: nEffFunc(wavelength)
        e = self.eps(wavelength)
        return e**0.5

        Material.__init__(self, wlMin, wlMax)

class FileWg(File):
    def __init__(self, filenameNEff, filenameBulk, lnNewBulkMaterial=None, mode=0, delimiter=','):
        txt1 = np.loadtxt(filenameNEff, delimiter=delimiter)
        wls = txt1[:,0] * 1.e-9
        nEffs = txt1[:,mode+1]
        wlMin = wls[0]
        wlMax = wls[-1]
        nEffFunc = spi.interp1d(wls, nEffs)
        self._eps = lambda wavelength: nEffFunc(wavelength)

        txt2 = np.loadtxt(filenameBulk, delimiter=delimiter)
        wls = txt2[:,0] * 1.e-9
        nBulks = txt2[:,1]
        nBulkFunc = spi.interp1d(wls, nBulks)

        nIndContrastFunc = lambda wls: nEffFunc(wls) - nBulkFunc(wls)
        if lnNewBulkMaterial:
            nEffNewFunc = lambda wls: lnNewBulkMaterial.n(wls*1.e9) + nIndContrastFunc(wls)
        else:
            nEffNewFunc = nEffFunc

        self.dn = lambda wls: nIndContrastFunc(wls*1.e-9)

        Material.__init__(self, wlMin, wlMax)

class Air(_Material):
    def __init__(self):
        _Material.__init__(self)

    def _eps(self, wavelength=None):
        return 1.

# http://www.opticsinfobase.org/view_article.cfm?gotourl=http%3A%2F%2Fwww.opticsinfobase.org%2FDirectPDFAccess%2FEEE59E78-F228-BD74-8B64990932FFCE71_69785%2Fao-41-24-5040.pdf%3Fda%3D1%26id%3D69785%26seq%3D0%26mobile%3Dno&org=Royal%20Melbourne%20Institute%20of%20Technology%20Swanston
class Ktp(_Material):
    def __init__(self, axis):
        _Material.__init__(self)
        assert(axis in ['x', 'y', 'z'])
        self.A = [None]*5
        if axis is 'x':
            self.A[0] = 3.29100
            self.A[1] = 0.04140
            self.A[2] = 0.03978
            self.A[3] = 9.35522
            self.A[4] = 31.45571
        elif axis is 'y':
            self.A[0] = 3.45018
            self.A[1] = 0.04341
            self.A[2] = 0.04597
            self.A[3] = 16.98825
            self.A[4] = 39.43799
        elif axis is 'z':
            self.A[0] = 4.59423
            self.A[1] = 0.06206
            self.A[2] = 0.04763
            self.A[3] = 110.80672
            self.A[4] = 86.12171

    # Permittivity
    def _eps(self, wavelength):
        A = self.A
        wavelengthUm = wavelength * 1.e-3 # [nm] -> [um]
        return A[0] + A[1] / (wavelengthUm**2 - A[2]) + A[3] / (wavelengthUm**2 - A[4])

# http://www.goochandhousego.com/wp-content/pdfs/LNmatProperties.pdf
class Ln(_Material):
    def __init__(self, axis, temperatureCelcius=20.):
        _Material.__init__(self)
        assert(axis in ['o', 'e'])
        self.T = temperatureCelcius
        self.F = (self.T - 24.5) * (self.T + 570.5)

        self.A = [None]*4
        self.B = [None]*3
        if axis is 'e':
            self.A[0] =  4.582
            self.A[1] =  9.921e4
            self.A[2] =  2.109e2
            self.A[3] =  2.194e-8
            self.B[0] =  5.2716e-2
            self.B[1] = -4.9143e-5
            self.B[2] =  2.2971e-7
        elif axis is 'o':
            self.A[0] =  4.9048
            self.A[1] =  1.1775e5
            self.A[2] =  2.1802e2
            self.A[3] =  2.7153e-8
            self.B[0] =  2.2314e-2
            self.B[1] = -2.9671e-5
            self.B[2] =  2.1429e-8

    # Permittivity
    def _eps(self, wavelength):
        A, B, F = self.A, self.B, self.F
        return A[0] + (A[1] + B[0]*F) / (wavelength**2 - (A[2] + B[1]*F)**2) + \
               B[2]*F - A[3]*wavelength**2

class Tfln(Ln):
    def __init__(self, axis, temperatureCelcius=20.):
        Ln.__init__(self, axis, temperatureCelcius)

        eps_1550_orig = super(Tfln, self)._eps(1550)

        if axis == 'o':
            no_1550 = 2.20600
            epso_1550 = no_1550**2
            self._deps = epso_1550 - eps_1550_orig
        elif axis == 'e':
            ne_1550 = 2.14455
            epse_1550 = ne_1550**2
            self._deps = epse_1550 - eps_1550_orig

    def _eps(self, wavelength):
        e = super(Tfln, self)._eps(wavelength)
        e += self._deps
        return e

class LnMg(_Material):
    def __init__(self, axis):
        _Material.__init__(self)
        assert(axis in ['o', 'e'])
        self.A = [None]*6
        if axis is 'e':
            self.A[0] = 2.2454
            self.A[1] = 0.01242
            self.A[2] = 1.3005
            self.A[3] = 0.05313
            self.A[4] = 6.8972
            self.A[5] = 331.33
        elif axis is 'o':
            self.A[0] = 2.4272
            self.A[1] = 0.01478
            self.A[2] = 1.4617
            self.A[3] = 0.05612
            self.A[4] = 9.6536
            self.A[5] = 371.216

    # Permittivity
    def _eps(self, wavelength):
        A = self.A
        wavelengthUm = wavelength * 1.e-3 # [nm] -> [um]
        e = wavelengthUm**2*(A[0]/(-A[1] + wavelengthUm**2)  + \
                             A[2]/(-A[3] + wavelengthUm**2)  + \
                             A[4]/(-A[5] + wavelengthUm**2)) + \
                             1.
        return e

# Gayer, 2008, Temperature and wavelength dependent refractive index equations for MgO-doped congruent and stoichiometric LiNbO3
class LnMgTemp(_Material):
    def __init__(self, axis, temperatureCelcius=20.):
        _Material.__init__(self)
        assert(axis in ['o', 'e'])
        self.T = temperatureCelcius
        self.F = (self.T - 24.5) * (self.T + 570.82)
        self.A = [None]*6
        self.B = [None]*6
        if axis is 'e':
            self.A[0] = 5.756
            self.A[1] = 0.0983
            self.A[2] = 0.2020
            self.A[3] = 189.32
            self.A[4] = 12.52
            self.A[5] = 1.32e-2
            self.B[0] = 2.860e-6
            self.B[1] = 4.700e-8
            self.B[2] = 6.113e-8
            self.B[3] = 1.516e-4
        elif axis is 'o':
            self.A[0] = 5.653
            self.A[1] = 0.1185
            self.A[2] = 0.2091
            self.A[3] = 89.61
            self.A[4] = 10.85
            self.A[5] = 1.97e-2
            self.B[0] = 7.941e-7
            self.B[1] = 3.3134e-8
            self.B[2] = -4.641e-9
            self.B[3] = -2.188e-6

    # Permittivity
    def _eps(self, wavelength):
        a = self.A
        b = self.B
        wavelengthUm = wavelength * 1.e-3 # [nm] -> [um]
        e = a[0] + b[0]*self.F + \
            (a[1] + b[1]*self.F) / (wavelengthUm**2 - (a[2] + b[2]*self.F)**2) + \
            (a[3] + b[3]*self.F) / (wavelengthUm**2 - a[4]**2) - \
            a[5]*wavelengthUm**2
        return e

# https://www.coherent.com/downloads/BBO_DS.pdf
class Bbo(_Material):
    def __init__(self, axis):
        _Material.__init__(self)
        assert(axis in ['o', 'e'])
        self.A = [None]*4
        if axis is 'e':
            self.A[0] = 2.3730
            self.A[1] = 0.0128
            self.A[2] = 0.0156
            self.A[3] = 0.0044
        elif axis is 'o':
            self.A[0] = 2.7405
            self.A[1] = 0.0184
            self.A[2] = 0.0179
            self.A[3] = 0.0155

    # Permittivity
    def _eps(self, wavelength):
        A = self.A
        wavelengthUm = wavelength * 1.e-3 # [nm] -> [um]
        e = A[0] + A[1]/(wavelengthUm**2 - A[2]) - A[3]*wavelengthUm**2
        return e

# https://books.google.com.au/books?id=zKI4hdtEVHwC&pg=PA216&lpg=PA216&dq=bibo+refractive+index+sellmeier&source=bl&ots=lyfLo24tVp&sig=z-shFbjI1HXynIKkS0XZENjoxOw&hl=en&sa=X&ei=KpNuVbjyJ8rz8gWk-YKIDg&ved=0CDoQ6AEwBQ#v=onepage&q=bibo%20refractive%20index%20sellmeier&f=false
class Bibo(Bbo):
    def __init__(self, axis):
        _Material.__init__(self)
        assert(axis in ['x', 'y', 'z'])
        self.A = [None]*5
        if axis is 'x':
            self.A[0] = 3.0722
            self.A[1] = 0.0324
            self.A[2] = 0.0315
            self.A[3] = 0.0133
        elif axis is 'y':
            self.A[0] = 3.1669
            self.A[1] = 0.0372
            self.A[2] = 0.0348
            self.A[3] = 0.0175
        elif axis is 'z':
            self.A[0] = 3.6525
            self.A[1] = 0.0511
            self.A[2] = 0.0370
            self.A[3] = 0.0226

# Linear optical characterization of chalcogenide glasses
# G. Boudebs, S. Cherukulappurath, M. Guignard, J. Troles, F. Smektala, F. Sanchez
class Chalcogenide(_Material):
    def __init__(self, chalcogenideType):
        _Material.__init__(self)
        self.chalcogenideType = chalcogenideType
        cauchyCoefs = { 'As2S3'        : [5.41, 0.20,  0.14] ,
                        'As2Se3'       : [7.56, 1.03,  0.12] ,
                        'GeSe4'        : [5.73, 0.80, -0.18] ,
                        'Ge10As10Se80' : [5.73, 0.80, -0.18] }
        self.A = cauchyCoefs[self.chalcogenideType]

    def _eps(self, wavelength):
        A = self.A
        wavelengthUm = wavelength * 1.e-3
        e = A[0] + A[1] / wavelengthUm**2 + A[2] / wavelengthUm**4
        return e

class SiO2(_Material):
    def __init__(self):
       _Material.__init__(self)

    def _eps(self, wavelength):
        x = wavelength * 1.e-3
        e = 1+0.6961663/(1-np.power(0.0684043/x,2))+0.4079426/(1-np.power(0.1162414/x,2))+\
            0.8974794/(1-np.power(9.896161/x,2))
        return e

class Su8(_Material):
    def __init__(self):
        _Material.__init__(self)
        self.coefs = [1.5525, 0.00629, 0.0004]

    def _eps(self, wavelength):
        wavelength /= 1000.
        n = _Material.cauchy_equation(wavelength, self.coefs)
        eps = n**2
        return eps

# https://refractiveindex.info/?shelf=main&book=Al2O3&page=Malitson-o
class Al2O3(_Material):
    def __init__(self, axis):
        _Material.__init__(self)
        assert(axis in ['o', 'e'])
        if axis == 'o':
            self._eps_e_o = self._eps_o
        elif axis == 'e':
            self._eps_e_o = self._eps_e

    def _eps_o(self, wavelength):
        x = wavelength * 1.e-3
        e = 1+1.4313493/(1-np.power(0.0726631/x,2))+ \
            0.65054713/(1-np.power(0.1193242/x,2))+ \
            5.3414021/(1-np.power(18.028251/x,2))
        return e

    def _eps_e(self, wavelength):
        x = wavelength * 1.e-3
        e = 1+1.5039759/(1-np.power(0.0740288/x,2))+ \
            0.55069141/(1-np.power(0.1216529/x,2))+ \
            6.5927379/(1-np.power(20.072248/x,2))
        return e

    def _eps(self, wavelength):
        return self._eps_e_o(wavelength)

# https://refractiveindex.info/?shelf=main&book=TiO2&page=Devore-o
class TiO2(_Material):
    def __init__(self, axis):
        _Material.__init__(self)
        assert(axis in ['o', 'e'])
        if axis == 'o':
            self._eps_e_o = self._eps_o
        elif axis == 'e':
            self._eps_e_o = self._eps_e

    def _eps_o(self, wavelength):
        x = wavelength * 1.e-3
        e = 5.913+0.2441/(np.power(x,2)-0.0803)
        return e

    def _eps_e(self, wavelength):
        x = wavelength * 1.e-3
        e = 7.197+0.3322/(np.power(x,2)-0.0843)
        return e

    def _eps(self, wavelength):
        return self._eps_e_o(wavelength)
