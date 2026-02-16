%% This script imports BrainProducts EEG files into EEGLAB and saves them as .set files for preprocessing
% Dependencies: EEGLAB on path with BrainProducts import plugin installed
% Set directories for EEGLAB and data files
% Last updated by Danielle Hewitt, 2nd Jan 2024

% Set path for EEGLAB
eeglab_path = 'C:\\Users\\michellel\\Desktop\\EEGLab\\eeglab2024.0'; %% change for the PC being used
addpath(genpath(eeglab_path));

% Specify subjects and main directory where data is stored
subjects = {'19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38'};
% subjects = {'20'}; %% make sure files are renamed from initials to numeric values in the form 'T1_[subjectID]';
mainDirectory = 'C:\\Data\\TouchStudy\\SA\\';

[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% Looping over all subjects
for i = 1:numel(subjects)
    currentSubject = subjects{i};

    % Looping over all EEG files
    currentDirectory = [mainDirectory 'T1_' currentSubject '\\'];
    EEGfiles =  dir(fullfile(currentDirectory, '*.EDF')); % finding all EEG files in subject folder
    for j = 1:length(EEGfiles)
        currentEEGFile = EEGfiles(j).name; % defining current EEG file for import

        % Getting name of current file to make setname
        setfiles = split(currentEEGFile, '.'); % removing the .vhdr
        setName = char(setfiles(1));
        setName = strrep(setName, '2023', '');

        % Load into EEGLAB and save as .set file
        EEG = pop_readedf([currentDirectory currentEEGFile]);
        EEG.setname = setName;
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0, 'setname', setName, 'gui', 'off');
        EEG = eeg_checkset(EEG);
        EEG = pop_saveset(EEG, 'filename', [setName '.set'], 'filepath', currentDirectory);
        [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

        disp(['Saved: ' currentEEGFile])
    end
end

disp('All done!'); %% if no files are generated, check the path is correct
