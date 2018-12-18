from . import _material_base as mb
from . import material as mat

class LnAni(mb._MaterialAni):
    def __init__(self, cut, temp=20):
        self._temp = temp
        self._cut = cut.lower()
        assert self._cut in ('x', 'z')

    @property
    def xx(self):
        if self._cut == 'x':
            axis = 'e'
        elif self._cut == 'z':
            axis = 'o'
        n = mat.Ln(axis, self._temp)
        return n

    @property
    def yy(self):
        if self._cut == 'x':
            axis = 'o'
        elif self._cut == 'z':
            axis = 'e'
        n = mat.Ln(axis, self._temp)
        return n

    @property
    def zz(self):
        if self._cut == 'x':
            axis = 'o'
        elif self._cut == 'z':
            axis = 'o'
        n = mat.Ln(axis, self._temp)
        return n
