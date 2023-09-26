# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 12:13:10 2023

@author: Lucian Spitzner


"""

import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm
import os
from os.path import exists 
import sys
import json
import resources.constants as const
import matplotlib.pyplot as plt
from joblib import Parallel, delayed

DIRPATH = os.path.dirname(os.path.abspath(__file__))

class plasmetheusSimulation:
    
    
    def __init__(self, setupFileName):
        """
    
        Parameters
        ----------
        setupFile : str
            name of setupFile in /data/setupFiles.

        Returns
        -------
        None.

        """
        
        self.simParams = {}
        self.fldParams = {}
        self.lineParams = {}
        self.simRes = {}
        
        with open(DIRPATH + '/setupFiles/' + setupFileName + '.json') as file:
            
            self.simParams = json.load(file)
            
            self.simParams['ns'] = len(self.simParams["specList"])


        # change certain values from SI to CGS
        
        self.simParams["stellarRad"] = self.simParams["stellarRad"] * const.R_sun


    def readLineList(self, species):
        
        def minFWHM(centre, deltaV=1):
            '''

            Parameters
            ----------
            centre : float
                centre of line in cm.
            deltaV : float, optional
                velocity dispersion in km/s. The default is 1.

            Returns
            -------
            float
                HWHM of absorption peak at centre with width of deltaV.

            '''
            
            # (c/lambda * delta_v/c = delta_v/lambda
            # multiply by 2 pi to get from FWHM to gamma 
            return 2*np.pi * (deltaV*1e5) / centre
        
        LineList = np.loadtxt(DIRPATH + '/resources/LineList.txt', dtype = str, usecols = (0, 1, 2, 3, 4), skiprows = 1)

        line_wavelength = np.array([x[1:-1] for x in LineList[:, 2]])
        line_A = np.array([x[1:-1] for x in LineList[:, 3]])
        line_f = np.array([x[1:-1] for x in LineList[:, 4]])
        
        SEL_COMPLETE = (line_wavelength != '') * (line_A != '') * (line_f != '') 

        SEL_SPECIES = (LineList[:, 0] == const.speciesInfoDict[species][0]) * (LineList[:, 1] == const.speciesInfoDict[species][1]) # Select element and ionization state 
        
        line_wavelength = line_wavelength[SEL_SPECIES * SEL_COMPLETE].astype(float) * 1e-8 # In cm
        original_gamma = line_A[SEL_SPECIES * SEL_COMPLETE].astype(float) #/ (4. * np.pi) # HWHM of Lorentzian, assumed to be A / (4 * pi) (see Draine 2013 p57)
        line_f = line_f[SEL_SPECIES * SEL_COMPLETE].astype(float)

        # see if the line width is wide enough (at least 1 km/s)
        min_gamma = minFWHM(line_wavelength)
        line_gamma = np.max((original_gamma, min_gamma), axis=0)
        
        SEL_WAVELENGTH = (line_wavelength > self.simParams["minWavelength"]) * (line_wavelength < self.simParams["maxWavelength"])

        n_incr = np.sum(np.array(original_gamma[SEL_WAVELENGTH]) < np.array(min_gamma[SEL_WAVELENGTH]))

        print(f'For {species} increased gamma in {n_incr} of {len(min_gamma[SEL_WAVELENGTH])} lines.')
        print('This is done to simulate spectraph resolutions.')
        
        self.lineParams[species] = [line_wavelength[SEL_WAVELENGTH], line_gamma[SEL_WAVELENGTH], line_f[SEL_WAVELENGTH]]

    
    def buildGrid(self):
        """
        builds the grid set by minWavelength and maxWavelength, as 
        well as resLow, resHigh and highResWidth. High-res regions are centered 
        in the original line position, and embedded in the uniform 
        low-res-region.
        
        Saves both the wavelength (wvl) and frequency (freq) grid.
        """
        # directly taken from Prometheus by Andrea Gebek 

        linesList = []
        
        for key_species in self.simParams["specList"]:
    
            # Read in center wavelengths of all absorption lines
            
            self.readLineList(key_species)
            
            linesList.extend(self.lineParams[key_species][0])
            

        peaks = np.sort(np.unique(linesList))
        diff = np.concatenate(([np.inf], np.diff(peaks), [np.inf])) # Difference between each absorption line center. Add infinity at beginning and end such that
        # when calculating the high-resolution bands we keep the lines with the lowest/highest wavlengths.
    
        HighResBorders = ([], [])
    
        for idx, peak in enumerate(peaks):
    
            if diff[idx] > self.simParams["highResWidth"]:
    
                HighResBorders[0].append(peak - self.simParams["highResWidth"] / 2.)
    
            if diff[idx + 1] > self.simParams["highResWidth"]:
    
                HighResBorders[1].append(peak + self.simParams["highResWidth"] / 2.)
    
        grid = []
        
        for idx in range(len(HighResBorders[0])):
    
            grid.append(np.arange(HighResBorders[0][idx], HighResBorders[1][idx], self.simParams["resHigh"]))
            
            if idx == 0 and self.simParams["minWavelength"] < HighResBorders[0][0]:
                grid.append(np.arange(self.simParams["minWavelength"], HighResBorders[0][0], self.simParams["resLow"]))
                
            if idx == len(HighResBorders[0]) - 1 and self.simParams["maxWavelength"] > HighResBorders[1][-1]:
                grid.append(np.arange(HighResBorders[1][-1], 
                                      self.simParams["maxWavelength"]+self.simParams["resLow"], 
                                      self.simParams["resLow"]))
                
                grid.append(np.arange(HighResBorders[1][idx - 1], HighResBorders[0][idx], self.simParams["resLow"]))
                
            else:
                grid.append(np.arange(HighResBorders[1][idx - 1], HighResBorders[0][idx], self.simParams["resLow"]))
    
        wavelength = np.sort(np.concatenate(grid))
        
        self.simParams["grid_wvl"] = wavelength
        self.simParams["grid_freq"] = const.c/wavelength
        self.simParams["gridLen"] = len(wavelength)
        

    def readFieldFile(self):
        """
        reads the AMITIS field file, adding the simulation box parameters to
        the dictionary self.fldParams
        
        Also reads the particle densities from the field file, turning them 
        into the required column densities. At this moment, only assumes
        single ionised species.


        Returns
        -------
        None.

        """
        
        def convStateNotation(speciesName):
            """
            

            Parameters
            ----------
            speciesName : str
                Elemental symbol, together with any number of "+", representing
                ionisation states. For example, single ionised oxygen is O+,
                double ionised magnesium Mg++

            Returns
            -------
            speciesName : str
                converted notation, where O+ becomes OII, Mg++ becomes MgIII.
                (i.e. OI is ground state).

            """
            
            numPlus = 0
            
            while speciesName.endswith('+'):
                
                speciesName = speciesName[:-1]
                
                numPlus += 1
                
            speciesName += 'I' + numPlus*'I'
            
            return speciesName
        
        with h5py.File(DIRPATH + '/data/' + self.simParams["dataFolder"] + "/" + self.simParams["fieldFileName"] + '.h5', 'r') as field:
            
            # List of all wanted attributes of the amitis field
            attributes = ['nx', 'ny', 'nz',
                          'dy', 'dz',       # in meter
                          'xmin', 'xmax',   # in meter
                          'ymin', 'ymax',   # in meter
                          'zmin', 'zmax',   # in meter
                          'obs_radius']     # in meter
            
            
            # add attributes to class dictionary
            for attr in attributes:
            
                try:
                    
                    self.fldParams[attr] = field.attrs[attr]
                    
                    # extract species names
                except KeyError:
                
                    raise KeyError(f'Field parameter {attr} could not be found in'+
                                   ' Amitis field attributes.')
                    
                    
            # add species ID used in AMITIS of all species 
            for keys in field.attrs.keys():
        
                if keys.startswith('s') and keys.endswith('_name'):
                    
                    speciesName = field.attrs[keys].decode('UTF-8')
                    
                    speciesName = convStateNotation(speciesName)
                    
                    # this removes 's' from the star and 
                    # the '_name' from the end.
                    # -1 to start counting from 0
                    speciesID = int(keys[:-5][1:]) - 1
                    
                    self.fldParams[speciesName] = speciesID
                    
        
        
            # read densities of selected species
            
            self.fldDens = {}
            
            for spec in self.simParams["specList"]:
            
                # add 1 since the field file count starts at 1, not at 0
                partID = self.fldParams[spec] + 1
                
                # read data and flatten the array 
                specChargeDens = np.array(field['rho0' + str(partID)]).flatten()
                
                # convert from 1/m**3 to 1/cm**3 and from charge density to number 
                # density and from density to column density
                # (assuming single ionized species)
                specDens = specChargeDens * 1e-6 / 1.602176634e-19 * (field.attrs['dx'] * 1e2)
                
                self.fldDens[spec] = specDens
            
        
    def createFilteredFile(self, hID=0):
        """
        Create a filtered particle file for faster file reading. It includes 
        three-dimensional position and velocity component in radial direction 
        for non-proton species. Removes hydrogen particles, as well as particles
        obstructed by the planet.

        Parameters
        ----------
        hID : int, optional
            ID of the SID of hydrogen in the AMITIS files. The default is 0.

        Returns
        -------
        None.

        """
    
    
        
        if exists(DIRPATH + "/data/" + self.simParams["dataFolder"] + "/" + self.simParams["partFileName"] + "_filtered.h5"):
            
            raise FileExistsError("Filtered AMITIS particle file already exists. Function call is likely bugged.")
                                  
        with h5py.File(DIRPATH + '/data/' + self.simParams["dataFolder"] + "/" + self.simParams["partFileName"] + '.h5') as partFile:
            
            # read SIDs in field file to exclude hydrogen
            sid = np.array(partFile['sid'])
            
            # Assuming that the sid of hydrogen is 0 in Amitis
            hydrogenMask = sid != hID
            
            newSid = sid[hydrogenMask]
            
            print('    Reading field dimensions')
            nx, ny, nz = self.fldParams['nx'], self.fldParams['ny'], self.fldParams['nz']
            xmin, xmax = self.fldParams["xmin"], self.fldParams["xmax"]
            ymin, ymax = self.fldParams["ymin"], self.fldParams["ymax"]
            zmin, zmax = self.fldParams["zmin"], self.fldParams["zmax"]
            pl_rad = self.fldParams["obs_radius"]
            
            
            # read positions of non-H particles
            rx = np.array(partFile['rx'])[hydrogenMask]
            ry = np.array(partFile['ry'])[hydrogenMask]
            rz = np.array(partFile['rz'])[hydrogenMask]
            
            # read x-components of velocities of non-H particles
            vx = np.array(partFile['vx'])[hydrogenMask]
            
            
            # construct the grid
            xrange, xstep = np.linspace(xmin, xmax, num=nx, endpoint=False, retstep=True)
            
            yrange, ystep = np.linspace(ymin, ymax, num=ny, endpoint=False, retstep=True)
            
            zrange, zstep = np.linspace(zmin, zmax, num=nz, endpoint=False, retstep=True)
            
            # find particles position in x,y,z index coordinates

            xpos = np.searchsorted(xrange, rx) - 1

            ypos = np.searchsorted(yrange, ry) - 1

            zpos = np.searchsorted(zrange, rz) - 1
            
            
            # flatten 3-dim index (voxels)
            loc_id = np.ravel_multi_index((xpos, ypos, zpos), (nx, ny, nz), mode='clip')
            
            # flatten 2-dim index (columns)
            loc_id_column = np.ravel_multi_index((ypos, zpos), (ny, nz), mode='clip')                         
            
            
            # exclude particles obstructed by planet
            vMask = np.sqrt(ry**2 + rz**2) > pl_rad
            
            print(f'Visible particles: {sum(vMask)}' + 
                  f'({sum(vMask)/len(rx)*100:.3f}% of all particles)')
            

            # create filtered h5py file
            
            with h5py.File(DIRPATH + '/data/' + self.simParams["dataFolder"] + "/" + self.simParams["partFileName"] + '_filtered.h5', 'w-') as newFile:

                newFile['sid'] = newSid[vMask]
                
                newFile['vx'] = vx[vMask]
                
                newFile['loc_id'] = loc_id[vMask]
                
                newFile['loc_id_column'] = loc_id_column[vMask]
    
        
    def readPartFile(self):
        """
        reads the filtered particle file, and will create a new one from the 
        unfiltered one if it doesnt exist.

        Parameters
        ----------
        folderName : str
            folder name containing desired AMITIS output
            in /data. Must contain a  field and (at least one) particle file.

        Returns
        -------
        None.

        """
        
        if not exists(DIRPATH + "/data/" + self.simParams["dataFolder"] + "/" + self.simParams["partFileName"] + "_filtered.h5"):
            
            self.createFilteredFile()

        with h5py.File(DIRPATH + '/data/' + self.simParams["dataFolder"] + "/" + self.simParams["partFileName"] + '_filtered.h5', 'r') as part:
            
            sid = np.array(part['sid'])
            loc_id = np.array(part['loc_id'])
            loc_id_column = np.array(part['loc_id_column'])
            vx = np.array(part['vx'])
            
            dat = pd.DataFrame(
                {
                    "vx": vx, 
                    "loc_id": loc_id,
                    "loc_id_column": loc_id_column,
                    "sid": sid
                }
            )
            
            self.partData = dat
    
    
    def setup(self):
        
        self.buildGrid()
        
        self.readFieldFile()
        
        self.readPartFile()
        
    def simulate(self):
        
        nparts = np.array([])
        
        all_tau = np.zeros((self.simParams['ns'], 
                            self.fldParams['ny'], 
                            self.fldParams['nz'], 
                            self.simParams['gridLen']))
        
        freq_x = self.simParams['grid_freq'][:, None, None]
                           
        for specID, spec in enumerate(self.simParams['specList']):
            
            print(f'Grouping by columns for species {spec}', end='...')
            
            # select species and group by column
            columnTable = (self.partData[self.partData['sid'] == self.fldParams[spec]]).groupby('loc_id_column')
            
            binwidths = np.zeros((columnTable.ngroups))

            print('done!')
            
            
            # getting line data and extending arrays for broadcasting
            centre_x = self.lineParams[spec][0][None, None, :]
            gamma_x =  self.lineParams[spec][1][None, None, :]
            eco_x   =  self.lineParams[spec][2][None, None, :]

            
            # calculating optical depth
            
            print('Calculating optical depth...')
            
            # parellelization:

            # list_of_CCFs, list_of_CCF_Es, list_of_T_sums = zip(*Parallel(n_jobs=NC,verbose=2*NT)(delayed(do_xcor)(i)))

            # define the function that will get parallelized: take one column as input, and calculate the optical depth per wavelength

            def calc_col_tau(idx, col_id, column):
                
                column_tau = np.zeros((self.simParams['gridLen']))
                                   
                for vid, voxel in column.groupby('loc_id'):
                    
                    if len(voxel['vx']) == 1:
                        
                        # convert velocity from m/s to cm/s
                        pv_m, v_m = np.array([1e-5]), voxel['vx'] * 10**2
                        
                        # arbitrarily setting bin width to 1 km/s (1e5cm/s)
                        deltaV = 1e5
                        
                        #(self.binwidths).append(deltaV)
                        
                        #nparts = np.append(nparts, 1)
                        
                       
                        vel_x = np.array(v_m)[None, :, None]
                        pv_x = pv_m[None, :, None]


                        tau_arr = deltaV * const.e**2 * np.pi / (const.m_e*const.c) * np.sum(
                                (
                                    (4*gamma_x) / 
                                    (16*np.pi**2 * (freq_x- (1-vel_x/const.c) * (const.c/centre_x))**2 + gamma_x**2) * pv_x*eco_x), 
                                axis=(2,1)
                                )
                        
                        
                        # multiply with the particle column density at the location (weight)
                        
                        column_tau += tau_arr * self.fldDens[spec][vid]
                    
                    else:
                        
                        # limit maximum bin width, also given in m/s
                        minv, maxv = np.min(voxel['vx']), np.max(voxel['vx'])
                        
                        nBinsMin = np.ceil((maxv - minv)/self.simParams['maxBinWidth'])
                        
                        nBins = np.max([self.simParams['velBins'], nBinsMin])
                    
                        
                        #setting number of bins, also converting m/s to cm/s
                        pv, v = np.histogram(np.float64(voxel['vx']*1e2), bins=int(nBins), density=True)
                    
                        
                        # ingore velocities that don't exist in the voxel
                        pv_mask = pv != 0

                        deltaV = (v[1] - v[0])
                        
                        #fix binwidth saving
                        # save binwidths for statistics 
                        #binwidths = np.append(binwidths, deltaV)
                        
                        # number of particles
                        #nparts = np.append(nparts, len(voxel['vx']))
                        
                        # velocities and their probabilities in the voxel
                        pv_m, v_m = pv[pv_mask], v[:-1][pv_mask] + deltaV/2

                        vel_x = v_m[None, :, None]
                        pv_x= pv_m[None, :, None]

                        # here the magic happens
                        tau_arr = deltaV * const.e**2 * np.pi / (const.m_e*const.c) * np.sum(
                                (
                                    (4*gamma_x) / 
                                    (16*np.pi**2 * (freq_x- (1-vel_x/const.c) * (const.c/centre_x))**2 + gamma_x**2) * pv_x*eco_x), 
                                axis=(2,1)
                                )
                                             
                        # multiply with the particle column density at the location (weight)
                        
                        column_tau += (tau_arr * self.fldDens[spec][vid])
                    
        
                
                (colY, colZ) = np.unravel_index(col_id, (self.fldParams['ny'], 
                                                         self.fldParams['nz']))
                
                # add the optical depth
                #all_tau[specID, colY, colZ] = column_tau

                return column_tau

        

        results = Parallel(n_jobs=50, verbose=1)(
            delayed(calc_col_tau)(idx,col_id,column) for idx, (col_id, column) in tqdm(enumerate(columnTable))
        )

        # missing multi-species implementation. save result in matrix "all_tau" for overlapping wavelengths
        # in different species (make it global?)
        print('Calculation done. Summing over all species and calculating absorption')
        
        #tot_tau = np.sum(all_tau, axis=0)
        
        absorption = np.e**(-np.stack(results, axis=0))

        # # mask empty columns
        # tot_mask = np.sum(np.abs(tot_tau), axis=2)
        
        # #number of non-empty columns
        # nCols = np.sum(tot_mask)
        
        # print(f'Number of non-empty columns: {nCols}')
        
        # surface area of a single column in cm**2
        columnArea = self.fldParams['dy'] * self.fldParams['dz'] * 1e4
        
        stellarArea = np.pi * self.simParams['stellarRad']**2
        
        planetArea = np.pi * (self.fldParams['obs_radius']*1e2)**2

        nCols = len(absorption)
        
        # relative area of planet, cloud and star
        tot_abs = ((np.sum(absorption, axis=0) * columnArea) +
               stellarArea - nCols*columnArea - planetArea)/ stellarArea 
        
        # save optical depth at highest absorbing wavelength
        # min, since the signal is lowest at that wavelength
        # maxabs = np.argmin(tot_abs)
        
        # maxtau = tot_tau[:,:,maxabs]
        
        print('Saving result')
        
        with h5py.File(DIRPATH + '/results/' + self.simParams["fieldFileName"] + '_res.h5', 'w-') as resFile:
            
            resFile['absorption'] = tot_abs
            
            #resFile['opticalDepth'] = maxtau
            
            resFile['wavelengths'] = self.simParams['grid_wvl']
            
            resFile['binwidths'] = binwidths
            
            resFile['partperVoxel'] = nparts

            resFile.attrs['velBins'] = self.simParams['velBins']
            
            resFile.attrs['stellarRadius'] = self.simParams['stellarRad']
            
            resFile.attrs['speciesList'] = self.simParams['specList']
        
        print('Plasmetheus exits now.')
        




        
#%%

mysim = plasmetheusSimulation('medium_dipole')
    
mysim.setup()

mysim.simulate()
