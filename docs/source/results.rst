Results
=======

The results will be saved in the results-directory. They are saved in the .h5-format. The following parameters are be accessible:

   'absorption':
      the realtive flux reduction at each wavelength specified in

   'wavelength':
      the grid of wavelengths

   'nParticles':
      statistics for the number of particles in each non-emtpy voxel

   'opticalDepth':
      absorption in every column at maximum absorbing wavelength

   

Optional Parameters
^^^^^^^^^^^^^^^^^^^

Depending on the parameters specified in the :ref:'setup-File'<setup>, the following fields are also available:

   'phaseAbs':
      absorption for every vertical slice of columns 

   'allAbs':
      absorption for every column (A lot of data!)