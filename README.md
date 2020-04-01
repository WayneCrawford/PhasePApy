# PhasePApy 

## Reference 

Chen,C., and A. A. Holland, (2016). PhasePApy: A Robust Pure Python Package for Automatic 
Identification of Seismic Phases, Seismological Research Letters, 87(6), DOI: 10.1785/0220160019.

			
## About

PhasePApy is a Pure Python program package developed by Chen Chen (c.chen@ou.edu; 
c.chen8684@gmail.com) under the guidance of Austin Holland.

Version 1.1: 2016	--	FBpicker, AICDpicker, KTpicker; 1D and 3D Associator

PhasePApy is a Seismic Phase Picker and Associator program package, which is 
designed to identify and associate picks of seismic microseismic events. 
Please reference the paper if you use this code.

Version 1.0 has been used by Oklahoma Geological Survey (OGS) for near real
time earthquake monitoring since 2014. Here, we only provide the version 1.1 to public,
because we think this one is more user-friendly and easy to manipulate. You might need 
to modify the code a bit for your requirements. Also, we still think the code has some 
space for further clean-up and improvement, so please feel free to modify and update it. 

If you have certain ideas to improve this package, please let us know.
If you are going to use this package for any commercial purpose, please let us know.
			

## Components 
For more information check out the wiki [https://github.com/austinholland/PhasePApy/wiki#welcome-to-the-phasepapy-wiki]

The PhasePApy package composes of two sub-packages: PhasePicker and Associator. These 
two sub-packages can work jointly or separately, depending on users’ requirements. 

### PhasePicker
The PhasePicker contains three pickers: FBpicker, AICDpicker, and KTpicker

### Associator
The Associator includes the 1D and 3D Associator modules which uses phase arrival look up
tables to associate earthquakes. The associator was tested for local to regional earthquake
 P- and S-phases, but the algorithm can be accommodated for any distance and possibly phase.
 
 Association is performed in the following steps:
   - Create a TravelTime table
   - Input picks to a Pick table
     - There's no way to tell it the phase, it figures it out in the following steps
   - Identify Candidates
     - assocND.LocalAssociator.id_candidate_events() identifies as candidates
       all of one stations Pick pairs that are offset by a value consistent with
       the range of S-P times in the traveltime table
     - Should maybe make another candidate builder that compares first arrivals
       across stations?
   - Associate Candidates
     - assocND.LocalAssociator.associate_candidates compares the projected origin
       times for all candidates: if there are more than `nsta_declare` stations
       with origin times separated by less than `assoc_ot_uncert`, then an Association
       is declared. I think there can be more than one Association (if there is
       another good cluster at a different time)
   - Add Single-Phase stations
     - `assocND.LocalAssociator.single_phase()` adds in stations without a good
       S-P time.  I haven't looked at this code yet.  I wonder if a final_stations
       criteria wouldn't help for Lucky Strike, which may not have many good S-P
       times but which often has a good set of P-phase 

## License
PhasePApy is released as public domain under the Creative Commons Public Domain Dedication [CC0 1.0] (https://creativecommons.org/publicdomain/zero/1.0/legalcode). 
The software was developed at the [Oklahoma Geological Survey](http://www.ogs.ou.edu).

## Examples 
There are examples in the test directory if one has the required libraries the examples
should run directly out of a git clone of PhasePApy.  

Because the associators require travel-time tables, we have included a data set and travel-
time tables in the git repository.  This is not ideal, but because PhasePApy does not 
include a travel-time calculator for velocity models this approach seemed like the best. 
For this reason the repository is larger than it needs to be.



## Required Libraries
Currently PhasePApy does not have a Python installer for module installation the files 
need to be placed in a proper location and the PYTHON_PATH or sys.path should point to where
the library exists.

The PhasePApy package relies on open libraries: 
+ Obspy 
+ Numpy
+ Scipy 
+ MatPlotLib (Optional) if you are going to plot results with plotting scripts in this 
package, you need matplotlib as well. Otherwise, the users can visualize it based on their own methods.	




