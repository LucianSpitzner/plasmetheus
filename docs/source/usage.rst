Usage
=====

Plasmetheus uses AMITIS files to produce absorption spectra of transiting exoplanetary exospheres.

.. admonition:: Disclaimer

      In order to use Plasmetheus, you have to have access to AMITIS files. These are provided by Shahab Fatemi, author of the code. 
      For copyright reasons and due to large file sizes, they're not included as a demo.

.. _installation:

Installation
------------

To use Plasmetheus, first downlaod it through git:

.. code-block:: console

   (.venv) $ git clone https://github.com/LucianSpitzner/plasmetheus.git


The important files and folders created are the following

plasmetheus.py
    the the main script,


setupFiles
    folder that contains your setupFiles. Multiple setupFiles for the same AMITIS simulation can exist.


data
    folder for subfolders for each AMITIS simulation. Note that subfolders must be created by the user.


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
For each run of Plasmetheus, you can configure the settings within the .json files.

outputName
   name of the plasmetheus simulation result file.

minWavelength and maxWavelength
   define the range of the simulation 

highResWidth
   defines the region width around each line with high resolution (line centre +/- highResWidth/2).

resLow and resHigh
   defines the resolution of the two regions. Usually, changing this is not necessary.


specList
   the list of species to analyze in the simulation. Be away that they have to match with species included in the AMITIS simulation.

stellarRad
   radius of the host star of the system

velBins
   minimum number of bins in velocity domain. No change necessary (more info at (add hyperlink to explanation))

maxBinWidth
   limits the bins width (and therefore sets a lower bound for the number of velocity bins per voxel). Unit is in m/s (defaults to
   10_000 m/s = 10 km/s)


dataFolder
   subfolder name of AMITIS simulation as specified by user

fieldFileName and partFileName
   names of the AMITIS simulation files (without the .h5 ending)

"nCores"
   number of cores to use. Due to overheading, a number larger than 30 cores leads to a slowdown and is not recommended.