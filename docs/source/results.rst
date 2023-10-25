Results
=======

The results will be saved in the results-directory. They are saved in the .h5-format. The following parameters are accessible:

   _`'absorption'`:
      the relative flux reduction at each wavelength

   'wavelength': in cm
      the grid of wavelengths

   'nParticles': array(int)
      statistics for the number of tracer particles in each non-emtpy voxel

   'binwidths': array(float)
      statistics for the histogram binwidths in each non-emtpy voxels

   _`'opticalDepth'`: 
      optical depth in every column at maximum absorbing wavelength

In `'absorption'`_, the result is processed and normalised to give the relative decrease of the stellar spectrum. 
In `'opticalDepth'`, the result is not yet calculated in reference to the star, but given in the dimensionless optical depth 
:math:`\tau`.

Optional Parameters
^^^^^^^^^^^^^^^^^^^

Depending on the parameters specified in the :ref:`setup-File <setup>`, the following fields are also available:

   'phaseAbs':
      absorption for every vertical slice of columns 

   'allTau':
      Optical depth for every column (A lot of data!)