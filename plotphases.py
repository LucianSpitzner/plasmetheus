#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import h5py
import sys
import plasmetheus
import os
from tqdm import tqdm

DIRPATH = os.path.dirname(os.path.abspath(__file__))

#%%

def plotPhases(plasSim):

    with h5py.File('results/' + plasSim.simParams['setupFileName'] + '_res.h5', 'r') as file:

        myres = np.array(file['phaseAbs'])
    
        wavel = np.array(file['wavelengths'])*1e8


    ny, nWvl = np.shape(myres)

    allSlices = np.zeros(shape=(ny*2-1,len(wavel)))
    
    # Area calculation
    # surface area of a single column in cm**2
    columnArea = plasSim.fldParams['dy'] * plasSim.fldParams['dz'] * 1e4
    
    # surface area of the star
    stellarArea = np.pi * plasSim.simParams['stellarRad']**2
    
    # ignore planet for now
    # surface area of the planet
    #planetArea = np.pi * (self.fldParams['obs_radius']*1e2)**2


    nz = plasSim.fldParams['nz']

    print('Calculating phase-wise absorption..')

    for y in tqdm(range(1,2*ny)):

        miny, maxy = np.max((0, y-ny)), np.min((y,ny))

        # nSlices here is the number of y-slices for each phase
        nSlices = 1 + maxy - miny
       
        selData = np.sum(myres[miny:maxy], axis=0)

        allSlices[y-1] = ((stellarArea - nSlices * nz * columnArea) + (columnArea * selData)) / (stellarArea)


    fig, ax = plt.subplots(figsize=(17,12))

    totmin = np.min(np.min(allSlices, axis=1)/np.max(allSlices, axis=1), axis=0)
    totmin = np.min(np.min(allSlices, axis=1)/np.max(allSlices, axis=1), axis=0)



    # Function to update the plot for each frame
    def update(frame):
        ax.clear()
        ax.plot(wavel, allSlices[frame]/np.max(allSlices[frame]))
        ax.set_ylim(totmin - 0.05*(1-totmin), 1 + 0.05*(1-totmin)) 
        
    # Create the animation
    num_frames = ny*2 - 1

    print("Rendering animation... (this can take a few minutes)")

    ani = animation.FuncAnimation(fig, update, frames=num_frames, repeat=False)

    # Save the animation as a GIF
    FFwriter = animation.writers['ffmpeg'](fps=30)
    
    ani.save(DIRPATH + '/figures/' + plasSim.simParams['setupFileName'] + '_phases.mp4', writer=FFwriter)



if __name__ == '__main__':
    
    if(len(sys.argv) != 2 ):

        raise ValueError("Arguments are not set correctly!")

    fileName = str(sys.argv[1]) 

    mysim = plasmetheus.plasSim(fileName)

    mysim.readFieldFile()

    plotPhases(mysim)

    


# %%
