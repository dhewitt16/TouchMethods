# ERP Analysis

Here we demonstrate a method to analyse event-related potentials following the offset of dynamic touch. The following preprocessing steps can be followed. In our example, this has been run after TF Analysis (though the instructions can be followed if running ERP analysis only).

## Importing the data into EEGLAB
1.	Run Import_Files_S2.m to get .set files if not already run

This step should have been completed during time-frequency (TF) analysis – no need to run again if so.  

2.	Run 1_EpochFiles_ERP.m

This epochs the files around the offset (-6 seconds pre-offset to +1 second postoffset), and does some simple preprocessing (rereference to common average, filter 0.1 to 30 Hz plus 48-52 Hz notch filter, downsample to 256 Hz).
ICA is not necessary here unless you want to run it again. The data is saved as one concatenated file for all conditions, for each participant.  
Change the directories for EEGLAB and the data directory to match those on the computer being used.  
If there are any issues with the script, run EEGLAB first.  

## Data cleaning
3.	Run 2_AddICAWeights.m  

This presumes that the ERD analysis has been run already. If not, or if different decomposition is required, uncomment to add the ICA step in the previous script.  
The script loads the set file from the ERD analysis after ICA decomposition. The files are loaded and the ICA weight matrix saved.  
The ICA weight matrix for each participant is loaded onto the ERP file from the previous step. This file is then saved for subsequent processing.  
All ERD set files must be stored in a directory within the current directory (cfg.dir) called SetFiles (i.e., not within subject folders – search for the name and copy them here first). Change the directories for EEGLAB and the data directory to match those on the computer being used.  

4.	Data cleaning

Load each data file individually. Reject independent components representing eye blinks and horizontal eye movements (should be the same components as in the ERD analysis). This should be done manually but can be checked using IClabel (settings >90% eye movement, as standard) – should generate the same results and can be scripted.  
Scroll through data for any noisy channels which should be interpolated – can use the channel spectra for confirmation of outliers.  
Save the data file as _ICA, incase the artefact rejection needs to be repeated.  
Use a method to reject noisy data. Here, we use the legacy EEGLAB rejection method to reject individual trials. Reject trials surpassing +/-125mv amplitude between -1 to +1 seconds (this length is based on analysis period – modify if different epochs are used). Reject trials with joint probability of over 7 / 7 standard deviations.  
Save the data file as _cleaned.  
These two steps can be scripted, particularly if the analysis needs to be repeated.  

## Analysis and Visualisation
5.	Run 3_ERPFigures.m.  

This does the main epoching and analysis steps.  The cleaned data is loaded and analysis (sub)epochs are computed using EEGLAB (-200 to 800ms pre-stimulus-offset). Baseline correction is carried out (-200 to 0s pre-stimulus-offset). The data for each epoch is saved so it does not need to be recomputed each time the script is run (after running the first time, change to ‘computeEpochs = 0’).  
Data is transformed from sensor data to current source density using FieldTrip. All data is saved as a big file so it does not need to be recomputed each time the script is run (after running the first time, change to ‘computeCSD = 0’).  
Figures are generated – grand average topoplots for each condition, grand average ERPs for each condition, and topoplots in selected time windows.  
If statsExports = 1, data is exported for statistical analysis. The data file will consist of 4 conditions x participant#, representing the mean amplitude during the plottingTimesSelected for the selectedElectrode/s. The results can be plotted using the BarPlotsTouch-Methods code in Jupyter Notebooks.  
Change paths for Fieldtrip, chanlocs file, currentDirectory. Ensure ‘wname’ matches the name from cleaned data in previous step. The baseline period can be amended, but should match the artefact rejection period from previous step.  
After computing CSD, change the ‘mostRecentFile’ to the name of the exported file – AllERP_csd_[date of export].mat.  
Electrodes for plotting and exports can be specified in selectedElectrode. If multiple electrodes are specified, these will be averaged.  

## Exporting the data for statistical analysis
6.	Run 4_BigExport.m  

The previous step generates one export per electrode/cluster and time point. To create one large file, run this, specifying the folderPath where the exports are located. The script will find all export files and combine them into one large table with subheadings of electrode number and time window based on the export file name.  

7.	Proceed to statistical analysis
