%This program will assemble data from 4 brushing conditions into one long file.
%The file will have triggers B1_1 and B2_1 (B1/2_2/3/4) to indicate the
%onsets and offesets of brushings in four different conditions
%By Andrej Stancak (2020), adapted by Danielle Hewitt (2023)

%Added: downsampling to 256 Hz to reduce the computational demands on ICA
%and to have neat 1-s intervals for time-frequency analysis.

%% This also removes practice trials due to noisy data in these sections - make sure to change if not 5

%% %Set paths for EEGLAB, codes AND FIELDTRIP
eeglab_path = '/Users/dhewitt/Analysis/eeglab2022.1/';
addpath(genpath(eeglab_path));

%Specify subjects and main directory where data is stored
cfg = [];
cfg.dir     = '/Users/dhewitt/Data/Touch/TouchStudy2/';
cfg.sub =  {'01', '02', '04', '05', '06', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17','18','19','20','21','22','23','24','25'};
%cfg.sub = {'01', '02'};
cfg.cond   = {'SLOW_OTHER', 'SLOW_SELF', 'FAST_OTHER', 'FAST_SELF'};

cfg.epoch = [-6 1]; % for touch offset ERPs

%==========================================================================

% Looping over all subjects
for iSub = 1:size(cfg.sub,2) %loop to run all subjects
    currentSubject = cfg.sub{iSub};
    currentDirectory = [cfg.dir 'T2_' currentSubject '/'];                   %!! please change the backslash to forwardslash

    for iCond = 1:4 %loop to run all conditions
        currentCond = cfg.cond{iCond};

        file2load = dir(fullfile(currentDirectory, char(['*' currentCond '.set'])));
        wname = [currentDirectory file2load.name];

        if exist(wname) == 0
            disp(['File ' wname ' does not exist']);
            return;
        end


        %find the event markers
        X = pop_loadset(wname);
        T1=0; %these are the beginnings and ends of blocks
        count = 0;

        for k=1:size(X.event,2)
            if strcmp('S  2',X.event(k).type) == 1 %% TRIAL   END
                count = count + 1;
                T1(count)=X.event(k).latency; %% there will be 40 plus the number of practice trials, unless these are removed
            end
        end

        X.event = '';
        for z=1:count
            X.event(z).latency  = T1(z);
            X.event(z).duration = 1;
            X.event(z).channel  = 0;
            X.event(z).type     = ['B' num2str(iCond)];
            X.event(z).code     = 'Stimulus';
            X.event(z).urevent  = z;
        end

        %start the epoching around event markers
        X.urevent = X.event;
        etype = {X.event(1).type};
        XE =  pop_epoch(X,etype,[cfg.epoch(1,1) cfg.epoch(1,2)]); %TRIAL DURATION

        if size(XE.epoch,2)>45 %indicating practice trials
            %deleting practice trials
            XE = pop_selectevent( XE, 'epoch',1:10 ,'select','inverse','deleteevents','off','deleteepochs','on','invertepochs','off');
            XE = eeg_checkset( XE );
        elseif size(XE.epoch,2)>40 %indicating practice trials
            %deleting practice trials
            XE = pop_selectevent( XE, 'epoch',1:5 ,'select','inverse','deleteevents','off','deleteepochs','on','invertepochs','off');
            XE = eeg_checkset( XE );
        else
            XE = eeg_checkset( XE );
        end

        eval(['EEG' num2str(iCond) ' = XE;']);

    end

    %==========================================================================

    %concatenate files
    MEEG = pop_mergeset(EEG1,EEG2,'keepall');
    MEEG = pop_mergeset(MEEG,EEG3,'keepall');
    MEEG = pop_mergeset(MEEG,EEG4,'keepall'); MEEG = eeg_checkset( MEEG );

    %doing some preproc
    MEEG = pop_reref( MEEG, []);  MEEG = eeg_checkset( MEEG ); %rereferencing to common av - should we insert ref el back in?
    MEEG = pop_eegfiltnew(MEEG, 'locutoff',0.1,'hicutoff',30); MEEG = eeg_checkset( MEEG ); %filter 0.1-30 hz
    MEEG = pop_eegfiltnew(MEEG, 'locutoff',48,'hicutoff',52,'revfilt',1);
    MEEG = pop_resample(MEEG,256);
    MEEG = pop_runica(MEEG,'runica');

    %==========================================================================

    pop_saveset(MEEG,'filename',['T2_' currentSubject '_allbrusherp_v4.set'],'filepath',currentDirectory);

    disp(['Subject T' currentSubject ' epoched, merged file saved after some preprocessing'])
end
