Confocal Polarized Monte Carlo (MATLAB) 
Polarized Monte Carlo model of a scanning confocal Mueller matrix microscope. Photons are launched from each scan point into a slab made of microspheres or birefringent media. Photons are detected through a geometric confocal gate: the objective NA and a pinhole back-projected to the focal plane. Running the four input polarization states (H, V, P (+45 degrees), R (right circular) gives the 16 signals needed to build a Mueller matrix at each focal depth. 

Created by Medical Photonics Laboratory, Florida International University

Contents 
1. Polarized_Confocal_MonteCarlo_Microspheres.m 
2. Polarized_Confocal_MonteCarlo_Birefringence.m
3. Helper functions (parsave*.m, rotSphi.m, updateU.m)

Requirements
  - MATLAB with the Parallel Computing Toolbox
  - MatScat (MATLAB File Exchange)

To start: 
1. Add this repository and dependencies to the Matlab path
2. Set the output folder: 'baseDir' (microspheres) or 'outputDir' (birefringence).  
3. For a first test, set Nphotons = 1e3 and run. Each depth folder ('D1,D2...' for microspheres and 'depth10um',... for birefringence) should then contain the Stokes files for each input state. 
4. Restore Nphotons and run the full simulation. 
