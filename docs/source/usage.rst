Usage
=====

Plasmetheus uses AMITIS files to produce absorption spectra of transiting exoplanetary exospheres.

.. admonition:: Disclaimer

      In order to use Plasmetheus, you have to have access to AMITIS files. These are provided by Shahab Fatemi, author of the code. 
      For copyright reasons and due to large file sizes, they're not included as a demo.

.. _installation:

Installation
------------

To use Plasmetheus, navigate to your desired location and download it through git:

.. code-block:: console

   $ git clone https://github.com/LucianSpitzner/plasmetheus.git


Setting up a virtual environment is recommended. For example, when using conda, do 

.. code-block:: console

   $ conda create --name plasenv
   $ conda activate plasenv

and install the requires packages

.. code-block:: console

   (plasenv) $ conda install numpy h5py pandas matplotlib tqdm joblib

.. note::
   Dependent on the package, try using :code:`(plasenv) $ pip install <packages>`.

The important files and directories are the following

plasmetheus.py
    the the main script


plotphases.py
   (provisional) script to plot phase-dependant absorption

setupFiles
    folder that contains your setupFiles. Multiple setupFiles for the same AMITIS simulation can exist

data
    folder for subfolders for each AMITIS simulation. Note that subfolders must be created by the user

figures
   save location for figures and animations.


.. _setup:

Setup
-----

AMITIS Simulation Data
^^^^^^^^^^^^^^^^^^^^^^

To use Plasmetheus, create a subdirectory in the data directory and move your AMITIS simulation files there. These consist of two
.h5 files: The field file (fld), which contains information about the magnetic field and the (charge) density in the simulation box, and
the particle file (par), which contains six-dimensional spatial and velocity information for tracer particles.
More information about the .h5 file format can be found at <https://docs.h5py.org/en/stable/>

SetupFile
^^^^^^^^^

For each run of Plasmetheus, you have to specify the simulation settings within a .json file in the setupFile-directory. 
For reference, there is a setupTemplate. It contains the following fields:

   outputName: str
      Name of the plasmetheus simulation result file

   minWavelength, maxWavelength: float
      in Angstrom (Å)
      Defines the range of the simulation. Within these limits, a grid is constructed with two different step-sizes: the
      general stepsize is `resLow`_, while in regions defined by `highResWidth`_ the stepsize will be `resHigh`_. The latter one must 
      be smaller, indicating a higher resolution.

   _`highResWidth`: float
      in Angstrom (Å)
      Defines the region width around each line with high resolution (line centre +/- highResWidth/2).

   _`resLow`, _`resHigh`: float
      in Angstrom (Å)
      Define the resolution of the two regions. Usually, changing this is not necessary.



   specList: list[str]
      The list of species to analyze in the simulation. Currently, Plasmetheus only works with singly ionised species. 
      The correct notation is <elemental symbol> + <ionisation state>, where neutral species get are I, singly ionised (i.e. C+) 
      are II, and so on. Internally, AMITIS uses a different notation, which is converted to the one used on Plasmetheus. 
      Be aware that the selected species have to match with species included in the AMITIS simulation. These can be found 
      inside the AMITIS field file.

   stellarRad: float
      in solar radii :math:`(1 R_{\odot} = 6.957e8\,m)`
      Radius of the host star of the system 

   semiMajorAxis: float
      in astronomical units (au)
      length of the semimajor axis of the planet-star system 
   orbitalPeriod
      in days (d)
      period of the orbit 



   velBins: int
      Minimum number of bins in the velocity domain. No change necessary (more info at (add hyperlink to explanation)).

   maxBinWidth: float
      in m/s (defaults to 10,000 m/s = 10 km/s)
      Limits the bins width (and therefore sets a lower bound for the number of velocity bins per voxel)

   customGamma: float
      value in km/s.
      increase of the intrinsic line width of transition to mimick both velocity distribution of the absorbers as well as telescope resolution
      


   dataFolder: str
      subdirectory of AMITIS simulation as inside data directory specified by user

   fieldFileName and partFileName: str
      names of the AMITIS simulation files (without the .h5 ending)

   nCores: int
      number of cores to use. Due to overheading, a number larger than 30 cores leads to a slowdown and is not recommended.

   savePhaseAbs: boolean
      boolean: if *true*, saves absorption for each column-slice (needed for phase-dependant plotting)

   saveCompleteTau: boolean
      boolean: if *true*, saves optical depth at every wavelength for all radial columns (Will cause a large result file size)

