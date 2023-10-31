#%%
import numpy as np
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import copy
import h5py
import sys
import plasmetheus
import os
from tqdm import tqdm
import time
import resources.constants as const

DIRPATH = os.path.dirname(os.path.abspath(__file__))

#%%

def plotPhases(plasSim):
    """
    Read the result of a Plasmetheus simulation and create a phase-dependant simulation of the spectrum.
    To run this, the "savePhaseAbs" parameter in the setupFile has to be set to true.

    """

    with h5py.File('results/' + plasSim.simParams['setupFileName'] + '_res.h5', 'r') as file:

        myres = np.array(file['phaseAbs'])
    
        wavel = np.array(file['wavelength'])*1e8

    # number of voxel across stellar disc:

    nStellar = int(np.floor(plasSim.simParams['stellarRad'] * 2 / (plasSim.fldParams['dy']*1e2)))

    stellarMid = np.ceil(nStellar/2)

    ny, nWvl = np.shape(myres)

    nPhases = nStellar + ny - 1


    # initialise result arrays for absorption and wavelength
    allSlices = np.zeros(shape=(nPhases,len(wavel)))

    allWV = np.zeros(shape=(nPhases, len(wavel)))
    

    # Area calculation
    # surface area of a single column in cm**2
    columnArea = plasSim.fldParams['dy'] * plasSim.fldParams['dz'] * 1e4
    
    # surface area of the star
    stellarArea = np.pi * plasSim.simParams['stellarRad']**2
    
    # ignore planet for now
    # surface area of the planet
    #planetArea = np.pi * (self.fldParams['obs_radius']*1e2)**2

    # index of midpoint (position of planet)
    planetMid = plasSim.fldParams['ymax']/plasSim.fldParams['dy']


    nz = plasSim.fldParams['nz']

    print('Calculating phase-wise absorption..', end=' ')

    for y in range(1, nStellar + ny):

        # index range of plasmetheus simulation result
        miny = np.max([0, ny - y])

        maxy = np.min([ny + nStellar - y, ny])

        # number of y-slices for each phase
        nSlices = maxy - miny
       
        selData = np.sum(myres[miny:maxy], axis=0)

        allSlices[y-1] = ((stellarArea - nSlices * nz * columnArea) + (columnArea * selData)) / (stellarArea)


        # dopplershift: change in the wavelength range due to orbital radial velocity

        shift = - (stellarMid + planetMid - y) * plasSim.fldParams['dy']*1e2
        
        # % of total orbit (assuming sin \approx tan for small phases < 5 deg)

        phaseangle = 2*np.pi * (shift / (2*np.pi* plasSim.simParams['semiMajorAxis']))
        
        orbitalVel = 2*np.pi*plasSim.simParams['semiMajorAxis'] / plasSim.simParams['orbitalPeriod']

        v_rad = np.sin(phaseangle) * orbitalVel
        
        allWV[y-1] = wavel*(1 + v_rad/const.c)

    print('done!')

    fig, ax = plt.subplots(figsize=(17,12))

    # stellar disc and cloud

    # stellar circle
    stellardisc = patches.Circle((0,0), radius=plasSim.simParams['stellarRad']/1e2, fill=False)

    cmap = plt.get_cmap('hot_r').copy()
    cmap.set_bad('white')

    # filter out zeros from density
    plasSim.fldDens['totDens'].T[plasSim.fldDens['totDens'].T == 0] = np.NaN

    totmin = np.min(np.min(allSlices, axis=(0,1)))
    totmax = np.max(np.max(allSlices, axis=(0,1)))


    # Function to update the plot for each frame
    def update(frame):

        # spectrum
        ax.clear()
        maxabs = allWV[frame][np.argmin(allSlices[frame])]  
        ax.axvline(maxabs, ls='--', alpha=.5, c='blue')
        for spec in plasSim.simParams['specList']:
            line = plasSim.lineParams[spec][0]
            ax.axvline(line*1e8, ls='--', alpha=.5, c='black')
        ax.plot(allWV[frame], allSlices[frame])
        ax.set_ylim(totmin - 0.05*(1-totmin), 1 + 0.05*(1-totmin)) 
        ax.set_xlim(np.min(wavel), np.max(wavel))
        
        # moving cloud
        axins = ax.inset_axes(
            [0.6,0.05, 0.38,0.5]
        )
        axins.add_patch(copy.copy(stellardisc))
        # set limits slightly larger than size of star
        axins.set_xlim((-plasSim.simParams['stellarRad']/1e2*1.1, plasSim.simParams['stellarRad']/1e2*1.1))
        axins.set_ylim((-plasSim.simParams['stellarRad']/1e2*1.1, plasSim.simParams['stellarRad']/1e2*1.1))

        shift = - (stellarMid + planetMid - frame) * plasSim.fldParams['dy']

        # plot density
        axins.imshow(plasSim.fldDens['totDens'].T,
                    origin='lower',
                    extent=[plasSim.fldParams['ymin']+shift,
                            plasSim.fldParams['ymax']+shift,
                            plasSim.fldParams['zmin'], 
                            plasSim.fldParams['zmax']],
                    vmin=0.01,
                    cmap=cmap)
        

    # Create the animation

    print("Rendering animation (this can take a few minutes)...", end='', flush=True)

    t1 = time.time()

    ani = animation.FuncAnimation(fig, update, frames=nPhases, repeat=False)

    # Save the animation as a GIF
    FFwriter = animation.writers['ffmpeg'](fps=30)

    ani.save(DIRPATH + '/figures/' + plasSim.simParams['setupFileName'] + '_phases.mp4', writer=FFwriter)

    t2 = time.time()

    print(f"done! Time elapsed: {t2-t1:.2f}")


if __name__ == '__main__':
    
    if(len(sys.argv) != 2 ):

        raise ValueError("Arguments are not set correctly!")

    fileName = str(sys.argv[1]) 

    mysim = plasmetheus.plasSim(fileName)

    mysim.readFieldFile()

    [mysim.readLineList(key_species) for key_species in mysim.simParams["specList"]]

    mysim.readFieldFile(saveTotDens=True)

    plotPhases(mysim)

    


# %%
