.. _running:

Running the code
----------------

Copy the setupFile and adjust the parameters as necessary. 
Make sure that you have the AMITIS-data available and in the correct directory.
To run the code, navigate into the plasmetheus folder using :code:`cd foo/plasmetheus`, and make sure your venv is activated.
Then, run the code with

.. code-block:: console

   (plasenv) $ python3 plasmetheus.py <setupFileName>

If the specific AMITIS simulation has not been analysed before, Plasmetheus will create a filtered particle file and place it into the
corresponding data folder. Then, it will start the calculations.
The code will update you about its progress. Usually, the most time-intensive part is the radiative transfer. The calculation time can vary between minutes 
and hours, but a progress bar will indicate the current status.

Results are placed into the results directory, named as specified in the setup file. 
