import audiolib.plotting as al_plt
import audiolib.tools as al_tls
import matplotlib.pyplot as plt
import numpy as np

from audiolib.elac import ElectroDynamic

fname = 'test_impedance.txt'
fname_added_mass = 'test_impedance_added_mass.txt'

# Txt has to be in mag-rad-format and the column-names in txt-file have to be:
# Freq	Z_abs	Z_rad
# (columns separated by \t)
z_dict = al_tls.read_intac_txt(fname)
z_dict_added_mass = al_tls.read_intac_txt(fname_added_mass)

added_mass = 0.1 # kg
Sd = 0.0363 # m²
c = 343 # m / s
rho = 1.204 # kg / m³

dayton_ls10 = ElectroDynamic(
    name = 'Dayton LS10-44',
    c=c,
    rho=rho,
    f_z = z_dict['Freq'],
    z = z_dict['Z_abs'],
    f_z_added_mass = z_dict_added_mass['Freq'],
    z_added_mass = z_dict_added_mass['Z_abs'],
    added_mass=added_mass,
    Sd = Sd,
)

dayton_ls10.plot_z_params()
dayton_ls10.print_ts()

