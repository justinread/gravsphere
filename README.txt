:::: GRAVSPHERE + BINULATOR ::::
Welcome to GravSphere and the binulator.

Use the binulator (binulator.py) to bin your discrete velocity 
and photometric data. It will then output all the files that 
GravSphere (gravsphere.py) needs to mass model these data.

For further details, please see the comments in the 
python source files.


:::: DEPENDENCIES ::::
To run these codes, you will need to install:
python3.x, numpy, scipy, matplotlib and emcee

The code has been verified to work with:
python 3.6.3
scipy version 0.19.1
numpy version 1.13.3
matplotlib version 2.1.0
emcee version 3.0rc2


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


:::: EXAMPLES ::::
Many examples are included for you to play with. You can run
all of the Milky Way "classical" dwarfs (and CVnI), two
Gaia Challenge mocks, and more. Please have a look in ./Data/
for these examples, and see:

binulator_initialise_<name>.py
gravsphere_initialise_<name>.py 

for example initialisation scripts for running these models.

To run an example: 

1. First set the Output directory, <output_base> in
   contsants.py

2. Next select the galaxy to analyse in binulator.py. For
   example to run the Gaia Challenge PlumCuspOm mock, 
   select: from binulator_initialise_PlumCuspOm import *

3. You will need to set up folders inside your output directory
   to store the output files. For PlumCuspOm, for example,
   this should be <out_base>/GCmock/PlumCuspOm/ for 1000 tracers.
   This folder structure is defined by <outfile> in:
   binulator_initalise_PlumCoreOm.py

4. You can now run python binulator.py

5. Once binulator has run, you will need to add further folders
   for the gravsphere output. Normally, you will want to run with
   "VirialShape" parameters, which requires:
   <out_base>/GCmock/PlumCuspOm/VirialShape/
   Propermotion output requires similarly:
   <out_base>/GCmock/PlumCuspOm/Propermotion/
   If using Propermotions *and* VirialShapes, then you will need:
   <out_base>/GCmock/PlumCuspOm/Propermotion/VirialShape/
   If using the "cosmo_cprior = 'yes'" option, this will be stored
   in a further subfolder, e.g.:
   <out_base>/GCmock/PlumCuspOm/VirialShape/CosmoC/

6. Next, select the same galaxy in gravsphere.py to run. For
   example to run the Gaia Challenge PlumCuspOm mock, 
   select: from gravsphere_initialise_PlumCuspOm import *

7. You can now run python gravsphere.py. You should run first with
   codemode = 'run'. Once the run has completed, you can run 
   again with codemode = 'plot' to produce some plots and output
   ASCII files. This output will appear in the output folder
   structure defined above.

8. Note that both binulator and gravsphere can take some time to run.
   On a fairly fast modern laptop, binulator can take up to an hour,
   while gravsphere will need to run overnight.



Justin Read | 01/02/21 

