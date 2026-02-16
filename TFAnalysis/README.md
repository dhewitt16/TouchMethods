EEG Analysis Pipeline

Brushing Data

Preparation and preprocessing  
1.	Get all data into one folder containing individual subject folders  
a.	Suggested: …/Data/TouchStudy/SA/  
b.	Rename all subject files in the form ‘T1_[1-number of subjects]’  
2.	Open MATLAB, start EEGLAB if not already in path  
3.	Import data into EEGLAB, save as a .set file – using ImportFilesSA.m  
a.	Then rename all .set files to include the block condition in the name based on the trials specified in .csv files from psychopy: [T1_[SUBJECTID]_SLOW_ARM.set, … SLOW_PALM, FAST_ARM, FAST_PALM]  

<img width="276" alt="Screenshot 2024-02-26 at 07 05 52" src="https://github.com/dhewitt16/TouchStudy/assets/122755414/96164de2-5b8b-4e1f-a252-5967660e3ea4">

5.	Epoch into trials and concatenate conditions into one file – using EpochFilesSA.m
a.	Loads file, finds the stimulus onset trigger (S  1) and splits into trials of -2.5 to +5.5 around that trigger
b.	Removes the practice trials (1st 5 trials in Study 1 – change to 10 trials for Study 2)
c.	Concatenate multiple blocks from same participant into one data set 
d.	Using combined dataset, rereference to common average
e.	Filter the data: high pass 1 Hz, low pass 70 Hz, notch 48-52 Hz (or equivalent 2Hz window either side of mains line noise) – check this is correctly done as the filtering has changed on EEGLAB
i.	Currently kept 50 Hz line noise assuming this is correct for SA – confirm with Sahba – can also use CleanLineNoise instead if there are any variations
f.	Save all as 1 combined .set file for each participant: T1_[SUBJECTID]_[SA]_allbrush.set
  
  
7.	Data cleaning – semi-automatic pipeline here  
a.	Eye blink removal – ICA
i.	Plot > Component activations (scroll) > set threshold to 40. Identify components showing blinks or horizontal eye movements (usually in the first 5).\
ii.	Can be checked with Plot > Component activations > in 2D > Components [1:10] and confirmed with Tools > Classify components using ICLabel > Label components
iii.	Tools > Remove components from data > specify components to remove (e.g., [1:3], [1], [1, 3])
b.	Open the data and look at it (Plot > Channel data (scroll) > set threshold to 65)
i.	 See if any electrodes need to be interpolated due to artefactual signal which affects large proportion of trials, e.g. if there is electrode pop or it is picking up lots of noise. Interpolate maximum of 10% of electrodes (6 with 64 channel system)

--- Save the data as T1_[SUBJECTID]_[CE/UK]_allbrush_ica.set ---
The next section runs data cleaning in parallel. A better method could be chosen but this is the one we have used so far.

c.	Data cleaning – mark trials showing large muscular or movement artefacts, particularly those where there is a difference between baseline and trial – maximum 25% trials rejected per block. Do this manually to get an idea of what is rejected, and then use semi-automated methods for consistency between participants and researchers.
i.	Legacy rejection method (File > Preferences > If set, show all menu items from previous EEGLAB versions > tick > Ok)
ii.	Tools > Reject data epochs > Reject using all methods:
1.	Find abnormal values: upper limit of 125 uV, lower limit -125 uV. Calc/Plot > Update Marks
2.	Find improbable data: Single channel limit 5 SD, All channel limit 5 SD. Calculate. If very high, some electrodes may need interpolating – check. Can also use limits of 7 SD for both limits – keep a record.
3.	Reject marked trials.
6.	Preprocessing complete – proceed to frequency analysis.

Frequency Analysis
1.	Using BrushSpa.m
a.	Convert from EEGLAB to Fieldtrip (make sure fieldtrip is in the path)
b.	Fieldtrip loads data into array of [participants x electrodes x frequency components x time x blocks]
c.	Fieldtrip gets all the trials for each block, and gets the power within each frequency band (within the foi specified – 1 to 70 Hz) within the trial for the specified time window of interest (-2.5 to 5.5 around stimulus onset)
d.	Relative power in the band calculated by specifying the baseline interval (2 to -1), and computing relative change versus the trial – baseline correction not carried out yet
e.	ERD is computed manually after power spectral decomposition, by subtracting baseline interval from the rest of the trial and multiplying by 100
f.	Relative power for all conditions saved in G struct, with separate structs for each condition – G.ERD[cond].powspctrm 63 (els) x 70 (freq) x 61 (time)
g.	Figures plotted if fplot.draw = 1
2.	Using BrushStat.m for UK or SA data
a.	Specify participants for analysis
b.	Select frequency band of interest – script needs to be run individually for each frequency band
c.	Modify the epochTime if necessary (cfg.epochTime) – currently 0.5 to 4.9 to avoid sharp changes at the start and end of trial
d.	Select if statistics should be computed. 1 = yes, 0 = no.
e.	Script will load all participant’s data and insert this into one D struct. This is saved as a .mat file when the script is ran, so next time it is ran (with the same number of participants) the data will be loaded from the file to save time
f.	Topographic maps (and a blank electrode map) are created for the specified frequency band and epoch time window.
g.	Time frequency figures are created for the specific electrode/s of interest
h.	If statistics are requested, a permutation analysis will be computed across all electrodes, for the specific frequency band and time window.
i.	Significant electrodes are shown in a topographic map.
ii.	Data exported to grand file with all electrodes for the specified frequency band and time window of interest - this file is used for further analysis
