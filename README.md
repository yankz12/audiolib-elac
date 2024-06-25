# Audiolib Electro-Acoustics Library

A small library to speed up programming with electroacoustic transducers

## Install

Be aware that this package requires the audiolib.plotting and audiolib.tools modules in order to work.
Install beforehand from here:

`https://github.com/yankz12/audiolib-tools`

and from here:

`https://github.com/yankz12/audiolib-plotting`

Now install elac: Download as zip and extract or clone the repository via

`git clone https://github.com/yankz12/audiolib-elac.git`

Install the package via pip installation into the directory

`pip install ./audiolib-elac`

Import the package e.g. via

`from audiolib.elac import ElectroDynamic`

and create objects like

`speaker = ElectroDynamic(...)`

## Application

The two practical applications are:

1.  After electrical input impedance measurement:

&nbsp;&nbsp; Give spectrum of electrical input impedance of the speaker without
and with added mass to the membrane (c, rho, f_z, z_abs, z_rad,
f_z_added_mass, z_added_mass, added_mass). The instance will
automatically calculate the Thiele-Small-parameters. The
TS-parameters can be extracted with ElectroDynamic.Qts,
ElectroDynamic.Mms, etc.
&nbsp;&nbsp; The modeled electrical input impedance spectrum, derived from the
TS-parameters, can be checked via f_es, z_es = ts_to_imp()

2.  Only TS-Parameters are known, no measurement present:

&nbsp;&nbsp; Give all known TS-parameters. The modeled electrical input
impedance, derived from the TS-parameters, can be checked via 
f_es, z_es = ts_to_imp()

Be aware that Lec is only modeled as a simple inductance, which will give 
wrong impedance info for high frequencies. Implementation of
semi-inductance-model yet to come: (http://www.cfuttrup.com/Thorborg_31.pdf)
