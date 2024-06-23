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
Lec_estimation_range = [500, 1000]

added_mass = 0.1 # kg
Sd = 0.0363 # m²
c = 343 # m / s
rho = 1.204 # kg / m³

dayton_ls10 = ElectroDynamic(
    name = 'Dayton LS10-44',
    c=c,
    rho=rho,
    f_z = z_dict['Freq'],
    z_abs = z_dict['Z_abs'],
    z_rad = z_dict['Z_rad'],
    f_z_added_mass = z_dict_added_mass['Freq'],
    z_added_mass = z_dict_added_mass['Z_abs'],
    added_mass=added_mass,
    Lec_estimation_range=Lec_estimation_range,
    Sd = Sd,
)

fig, ax = dayton_ls10.plot_z_params()
dayton_ls10.print_ts()

f_es, z_es = dayton_ls10.ts_to_imp()
al_plt.plot_rfft_freq(
    f_es,
    np.abs(z_es),
    xscale='log',
    fig=fig,
    ax=ax,
    color='pink',
    linestyle='--',
)
