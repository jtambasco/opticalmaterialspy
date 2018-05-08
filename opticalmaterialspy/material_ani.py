from . import _material_base as mb
from . import material as mat

class Ln(mb._MaterialAni):
    def __init__(self, temp=20):
        self._temp = temp

    @property
    def xx(self):
        return mat.Ln('e', self._temp)

    @property
    def yy(self):
        return mat.Ln('e', self._temp)

    @property
    def zz(self):
        return mat.Ln('o', self._temp)
