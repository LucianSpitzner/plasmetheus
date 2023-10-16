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
   Dependent on the package, try using :code:`(plasenv) $ git install`.

The important files and folders are the following

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

Data
^^^^

To use Plasmetheus, create a subdirectory in the data folder and move your AMITIS simulation files there. These consist of two
.h5 files: The field file (fld), which contains information about the magnetic field and the (charge) density in the simulation box, and the
particle file (par), which contains six-dimensional spatial and velocity information for tracer particles. 

SetupFile
^^^^^^^^^
For each run of Plasmetheus, you can configure the settings within the .json files:

   outputName
      Name of the plasmetheus simulation result file

   minWavelength, maxWavelength
      Define the range of the simulation.

   highResWidth
      Defines the region width around each line with high resolution (line centre +/- highResWidth/2).

   resLow, resHigh
      Define the resolution of the two regions. Usually, changing this is not necessary.



   specList
      The list of species to analyze in the simulation. Be aware that they have to match with species included in the AMITIS simulation.

   stellarRad
      Radius of the host star of the system in solar radii :math:`(1 R_{\odot} = \num{6.957e8} \mathrm{ m})`

   semiMajorAxis
      length of the semimajor axis of the planet-star system in astronomical units (au)

   orbitalPeriod
      period of the orbit in days (d)



   velBins
      Minimum number of bins in the velocity domain. No change necessary (more info at (add hyperlink to explanation)).

   maxBinWidth
      Limits the bins width (and therefore sets a lower bound for the number of velocity bins per voxel). Unit is in m/s (defaults to 10,000 m/s = 10 km/s).

   customGamma
      increase of the intrinsic line width of transition to mimick both velocity distribution of the absorbers as well as telescope resolution
      value in km/s. best kept at 1


   dataFolder
      subfolder name of AMITIS simulation as specified by user

   fieldFileName and partFileName
      names of the AMITIS simulation files (without the .h5 ending)

   nCores
      number of cores to use. Due to overheading, a number larger than 30 cores leads to a slowdown and is not recommended.

   savePhaseAbs
      boolean: if true, saves absorption for each column-slice (needed for phase-dependant plotting)

   saveCompleteAbs
      boolean: if true, saves absorption for every radial column (:note:`Will cause a large result file size`)





.. _running:

Running the code
----------------

Copy the setupFile and adjust the parameters as necessary. 
To run the code, navigate into the plasmetheus folder using :code:`cd foo/plasmetheus`.
Then, run the code with::

   python3 plasmetheus.py <setupFileName>

If the specific AMITIS simulation has not been analysed before, Plasmetheus will create a filtered particle file and place it into the
corresponding data folder. Then, it will start the calculations.

Results are placed into the results directory, named as specified in the setup file. 

