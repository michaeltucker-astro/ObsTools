# ObsTools

Strongly suggest creating a new Anaconda environment for these tools, will prevent any dependency conflicts
Assuming you have Anaconda installed, do the following:
(replace ObsTools with whatever you want to name this environment)
$conda create -n ObsTools python=3.6 numpy scipy matplotlib astropy argparse pip
Then activate the environment:
$source activate ObsTools (or whatever you named it)
Now you're ready to roll.

If you feel a new conda env is unnecessary, then the same packages/python version are needed. Otherwise, might crash, might not, continue at your own peril.

To setup, run:
  $python setup.py
 
#will install any other necessary packages, and will ask if you want to make the files executable and/or add to $PATH
After that, should be good to go!

Package contents:

--extractSNIFS.py: extracts and combines SNIFS blue and red spectra into ASCII format, with comment headers for MJD-OBS, EXPTIME, etc.
--quick-class.py: tools for eyeballing redshift and line ID-ing of transient spectra. comments are encouraged
--obstools.py: simple observing tools, such as finder chart and converting CSV file to a file readable by JSkyCalc

To-do list:
--SNIFSphot.py: input SNIFS acq image, determine WCS solution, and apply aperture photometry using Pan-STARRS data
--obstools.py -> scheduler: given a target list, find optimal observing schedule to minimize telescope movement/delays
--transient-spec.c: C code to cross-corr observed spectra with templates, basically SNID with up-to-date templates
  -> will eventually implement machine-learning (required for fellowship and thesis)

If any suggestions/bugs/ideas/comments/improvements/anything else, let me know!
