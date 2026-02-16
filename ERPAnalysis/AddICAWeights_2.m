%% Touch Study 2 - ERPs - Adding ICA Weights (Step 3) %%

% This loads the ICA weight matrix generated for each participants' TF data.
% Then, it loads each subjects' data following epoching for ERPs, and adds
% the ICA weight matrix for further preprocessing. Next step: artefact
% rejection and channel interpolation.

% Danielle Hewitt (07/11/2023)

%% Set paths for EEGLAB
eeglab_path = '/Users/dhewitt/Analysis/eeglab2021.0/';
addpath(genpath(eeglab_path));

%==========================================================================

%Specify subjects and main directory where data is stored
cfg = [];
cfg.dir     = '/Users/dhewitt/Data/Touch/TouchStudy2/';
%cfg.sub =  {'01', '02', '04', '05', '06', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17','18','19','20','21','22','23','24','25'};
cfg.sub =  {'01'};
%==========================================================================

%% Loading the data for all subjects

for iSub = 1:size(cfg.sub,2)
    currentSubject = cfg.sub{iSub};
    currentDirectory = [cfg.dir 'SetFiles/'];

    wname = ([currentDirectory 'T2_' currentSubject '_allbrush.set']);

    if exist(wname) == 0
        disp(['File ' wname ' does not exist']);
        return
    end

    EEG = pop_loadset(wname);
    EEG = eeg_checkset( EEG );

    %==========================================================================

    currentDirectory = [cfg.dir 'T2_' currentSubject '/'];
    ica_wname = ([currentDirectory '/T2_' currentSubject '_combinedICAmatrix']);
    pop_expica(EEG, 'weights', ica_wname);

    %==========================================================================

    %% Next, applying that to each participant

    wname = ([currentDirectory 'T2_' currentSubject '_allbrusherp_v4.set']);

    if exist(wname) == 0
        disp(['File ' wname ' does not exist']);
        return
    end

    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);
    eeglab redraw;
    EEG = pop_loadset(wname);
    EEG = eeg_checkset( EEG );

    %epoching
    EEG = pop_editset(EEG, 'run', [], 'icaweights', ica_wname, 'icasphere', []);
    EEG = eeg_checkset( EEG );
    pop_saveset(EEG,'filename',['T2_' currentSubject '_allbrusherp_icaweights_v4.set'],'filepath',currentDirectory); %saving the file ready for ICA rejection and artifact rejection

    disp(['Subject T' currentSubject ' epoched, merged file saved after some preprocessing'])

end


