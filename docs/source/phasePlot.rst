
.. _phaseplot:

phase-dependant animations
--------------------------

.. note::
   Preliminary functionality:

A phase-dependant animation of the result can be created using :code:`python3 plotphases.py <setupFile>`.
For this, the option "savePhaseAbs" must be set to true. It will create an animation in the figure-dir.
This will calculate the relative flux reduction as the planet passes infront of the host star. For this, several approximation are made:

   - The planet does not move on a circular arc, but on a straight line instead. This approximation is more accurate for planets with a large semimajor axis.

   - The exopsheric cloud is not rotated. This will have to be implemented. Common rotations are for planets during transit are up to [-10,10] degrees.

   - The planetary shadow is currently ignored. 