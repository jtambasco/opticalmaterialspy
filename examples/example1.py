import opticalmaterialspy as mat

m = mat.SiO2()

# Refractive index @ 1550nm.
print('n(1.55e-6m):', m.n(1.55e-6)) # Knows 1.55e-6 must be [m].
print('n(1.55um):', m.n(1.55)) # Knows 1.55 must be [um].
print('n(1550nm):', m.n(1550)) # Knows 1550 must be [nm].

# Group velocity refractive index @ 900nm.
print('n_gv(900nm):', m.ng(900))

# Group velocity dispersion @ 808nm.
print('GVD(0.808um):', m.gvd(0.808))
