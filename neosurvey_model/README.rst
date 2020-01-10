NEOSurvey Model
===============

This is a backup of the original code that was used in the modeling of the
SpitzerNEOs targets. The model is described in detail by
`Trilling et al. (2016)
<https://ui.adsabs.harvard.edu/abs/2016AJ....152..172T/abstract>`_.

Note that this code is not cleaned in any way and looks like a mess. It is only
provided here for reference.

Compile
-------

>>> gcc neosurvey_model.c -o neosurvey_model -lm

Use
---

The following configuration has been used on all SpitzerNEOs data to derive
diameter and albedo estimates, including corresponding uncertainties:

>>> neosurvey_model example_data.dat 4n -montecarlo

You can run the model without error estimation using:

>>> neosurvey_model example_data.dat 4n

resulting in the following output::

    2201_Oljato_(1947_XC) (NEATM, fixed eta): d=2.468, pv=0.156, eta=1.266, chi2=0.000e+00
    137108_(1999_AN10) (NEATM, fixed eta): d=0.527, pv=0.316, eta=1.238, chi2=0.000e+00
    506425_(2000_DQ110) (NEATM, fixed eta): d=4.905, pv=0.013, eta=1.372, chi2=0.000e+00
    (2000_RN12) (NEATM, fixed eta): d=0.735, pv=0.030, eta=1.363, chi2=0.000e+00
    452376_(2002_AC5) (NEATM, fixed eta): d=0.244, pv=0.250, eta=1.692, chi2=0.000e+00
    (2002_LW) (NEATM, fixed eta): d=0.079, pv=0.216, eta=1.670, chi2=0.000e+00
    99942_Apophis_(2004_MN4) (NEATM, fixed eta): d=0.384, pv=0.245, eta=1.528, chi2=0.000e+00
    340048_(2005_VT5) (NEATM, fixed eta): d=0.328, pv=0.282, eta=1.332, chi2=0.000e+00
    378358_(2007_LD) (NEATM, fixed eta): d=0.314, pv=0.336, eta=1.292, chi2=0.000e+00
    (2007_TL23) (NEATM, fixed eta): d=0.135, pv=0.201, eta=1.448, chi2=0.000e+00
    (2010_WQ7) (NEATM, fixed eta): d=0.802, pv=0.104, eta=1.315, chi2=0.000e+00


Keep in mind that the results of the error estimation are stochastic, so
your results will deviate minimally from published results.


Input Data
----------

The input data file has to contain the following information in
this order
and a tabular way, separated by whitespaces. Each line corresponds to an
individual single-band observation.

* target identifier (string)
* distance from the observer (float, au)
* heliocentric distance (float, au)
* solar phase angle (float, deg)
* Solar System absolute magnitude (float, mag, V-band)
* uncertainty of the absolute magnitude (float, mag)
* photometric slope parameter (float, H-G system)
* thermal observation wavelength (float, micron)
* thermal observation flux density (float, mjy)
* thermal observation flux density uncertainty (float, mjy)

License
-------

All code in this repository is distributed under a 3-clause BSD licence.