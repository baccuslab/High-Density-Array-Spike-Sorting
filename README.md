# High-Density-Array-Spike-Sorting
Spike sorting code for the high density array, provided by Felix Franke.

## Description of spike-sorting-tests

In the folder, spike-sorting-tests, the example script in the subdirectory Mysort/+mysort/+sorters/sortTest_Aran.m, 
uses a simulated data file (not included in this repo) from Espen Hagen (from the Gaute Einevoll lab) to test the spike 
sorting code. It is a slight modification of sortTest.m to read the data file locally. The data set from Espen is a simulation of a laminar electrode 
with 16 single electrodes. They are close enough that you can see the same neuron on 3-4 electrodes. 
As it is meant to be a simple test scenario, this script tests only the actual sorter for a small number of electronde.
In fact, only 8 of the electrodes into the sorter in this test, so, roughly, only half of the neurons are "seen" by the sorter.
The matrix you will get as an output are the sorting errors. As one can see, some neurons are very nicely sorted, some aren't.

## Description of spike-sorting-code

In the folder, spike-sorting-code, are the scripts used to sort the high density data. It uses the sorter used in spike-sorting-tests, 
however, to sort individual electrode groups and then merges, in a second step, all the individual electrode group sortings into one final 
consistent sorting. As a disclaimer, it is a bit difficult to use now, however, because the most recent incarnation of it relies on Felix's group's 
computation grid.

To start the sorting, Felix uses for each person (different ways of organizing their data) a different script. These scripts are in 
spike-sorting-code/DataAnalysis/+ana/+sorter/dataScripts

Only a small number of them are up to date. One script that is currently in use is:
startHillerMany2013_grid.m

These scripts then call another script in the same folder, entitled, gridSpikeSorting.m.

This in turn handles the preprocessing of the data and then submits the computation to the grid. The whole grid business is handled by this class:
spike-sorting-code/GridComputing/+grid/sortJob.m,

with the exception of the actual submit process which Felix implemented simply via the file system 
(there is a matlab process checking a folder for files. If it finds a file, it moves it somewhere else and uses this file to look in the data folder for a prepared grid job).

The next entry point is then on each of the grid nodes in the SortJob class. In the function "run", the actual sorting of an 
electrode group is started. The sortJob is the class that handles the splitting up into tasks for their grid engine.
"run" is then executed for each task. In the function "postprocessBOTMSorting", after all tasks are done by the grid, the results 
are collected at a single location and merged into a final sorting.



