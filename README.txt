:::: GRAVSPHERE + BINULATOR ::::
Welcome to GravSphere and the binulator.

Use the binulator (binulator.py) to bin your discrete velocity 
and photometric data. It will then output all the files that 
GravSphere (gravsphere.py) needs to mass model these data.

For further details, please see the comments in the 
python source files.


:::: DEPENDENCIES ::::
To run these codes, you will need to install:
python3.0, numpy, scipy, matplotlib and emcee


:::: CITATIONS ::::
If using this code, please cite the following code papers:

https://ui.adsabs.harvard.edu/abs/2017MNRAS.471.4541R/abstract
https://ui.adsabs.harvard.edu/abs/2018MNRAS.481..860R/abstract
https://ui.adsabs.harvard.edu/abs/2020MNRAS.498..144G/abstract

If using the J-factor integral calculation, please cite also:
https://ui.adsabs.harvard.edu/abs/2020JCAP...09..004A/abstract

Please acknowledge use of the gravsphere code and/or
binulator, and link to its github page:
https://github.com/justinread/gravsphere


:::: EXAMPLES ::::
Many examples are included for you to play with. You can run
all of the Milky Way "classical" dwarfs (and CVnI), two
Gaia Challenge mocks, and more. Please have a look in ./Data/
for these examples, and see:

binulator_initialise_<name>.py
gravsphere_initialise_<name>.py 

for example initialisation scripts for running these models.


:::: BUGS ::::
If you spot any bugs, please let me know!


:::: NOTES ::::
Note that this public release of the code is distinct from
the "FreeForm" version used in most of the above publications. 
In particular, it bins the data in an entirely new way that is
more robust for smaller datasets and "ultra-faints", where 
the velocity errors approach the velocity dispersion. This 
yields slightly larger uncertainties, however, than the method
presented in the original GravSphere code paper:
https://ui.adsabs.harvard.edu/abs/2017MNRAS.471.4541R/abstract

If you would like to use the fully free form version, with
the original binning methodology, please see Anna Genina's 
independent public release here: 
https://github.com/AnnaGenina/pyGravSphere


Justin Read

