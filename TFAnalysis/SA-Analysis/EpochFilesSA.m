%This program will assemble data from 4 brushing conditions into one long file.
%The file will have triggers B1_1 and B2_1 (B1/2_2/3/4) to indicate the
%onsets and offesets of brushings in four different conditions
%By Andrej Stancak (2020), adapted by Danielle Hewitt (2023)
% Last updated: Danielle Hewitt, 3rd Jan 2024

%Added: downsampling to 256 Hz to reduce the computational demands on ICA
%and to have neat 1-s intervals for time-frequency analysis.

%make sure the data is renamed to include condition blocks in the setfile

%% This also removes practice trials due to noisy data in these sections - make sure to change if not 5

%% %Set paths for EEGLAB, codes AND FIELDTRIP
eeglab_path = '/Users/dhewitt/Analysis/eeglab2022.1/';
addpath(genpath(eeglab_path));

%Specify subjects and main directory where data is stored
cfg = [];
cfg.dir     = '/Users/dhewitt/OneDrive - Nexus365/Data/Touch/SAData24/'; % change to backslash for windows os if necessary
%cfg.sub  = {'05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15'};
cfg.sub = {'19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38'}
cfg.cond   = {'SLOW_PALM', 'SLOW_ARM', 'FAST_PALM', 'FAST_ARM'};

cfg.epoch = [-2.5 5.5]; %this epoch will suit a baseline such as -2 to -1 and to explore the post-brush rebound of beta components. Basically, the effective epoch is -2 to +5s


%==========================================================================

% Looping over all subjects
for iSub = 1:size(cfg.sub,2) %loop to run all subjects
    currentSubject = cfg.sub{iSub};
    currentDirectory = [cfg.dir 'T1_' currentSubject '/'];                   %!! change to backslash for windows os

    for iCond = 1:4 %loop to run all conditions
        currentCond = cfg.cond{iCond};

        file2load = dir(fullfile(currentDirectory, char(['*' currentCond '.set'])));
        wname = [currentDirectory file2load.name];

        if exist(wname) == 0
            disp(['File ' wname ' does not exist']);
            return;
        end

        % Renaming the EEG triggers
        % Firstly, find the event markers
        X = pop_loadset(wname);

        for i = 1:length(X.event)
            % Check if the 'type' column contains the specified string
            if strcmp(X.event(i).type, 'KB-Marker-m (S1) ') %trial onset
                % Replace the string with 'S 1'
                X.event(i).type = 'S  1';
            end
        end


        T1=0; %these are the beginnings and ends of blocks
        count = 0;

        for k=1:size(X.event,2)

            if strcmp('S  1',X.event(k).type) == 1 %% TRIAL START
                count = count + 1;
                T1(count)=X.event(k).latency; %% there will be 40 plus the number of practice trials, unless these are removed
            end
        end

        % for k = 1:length(X.event)
        %     if strcmp('''KB-Marker-m (S1)''', X.event(k).type) == 1 %% TRIAL START
        %         count = count + 1;
        %         T1(count) = X.event(k).latency;
        %     end
        % end

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

        %deleting practice trials
        XE = pop_selectevent( XE, 'epoch',[1:5] ,'select','inverse','deleteevents','off','deleteepochs','on','invertepochs','off');
        XE = eeg_checkset( XE );

        eval(['EEG' num2str(iCond) ' = XE;']);

    end

    %==========================================================================

    %concatenate files
    MEEG = pop_mergeset(EEG1,EEG2,'keepall');
    MEEG = pop_mergeset(MEEG,EEG3,'keepall');
    MEEG = pop_mergeset(MEEG,EEG4,'keepall'); MEEG = eeg_checkset( MEEG );

    %doing some preproc
    MEEG = pop_editset(MEEG, 'chanlocs', [pwd '/gtec64.sph']);
    MEEG=pop_chanedit(MEEG, 'nosedir','+Y','rplurchanloc',1);
    MEEG = pop_reref( MEEG, []);  MEEG = eeg_checkset( MEEG ); %rereferencing to common av - should we insert ref el back in?
    MEEG = pop_eegfiltnew(MEEG, 'locutoff',0.1,'hicutoff',70); MEEG = eeg_checkset( MEEG ); %filter 1-70 hz
    MEEG = pop_eegfiltnew(MEEG, 'locutoff',48,'hicutoff',52,'revfilt',1);  MEEG = eeg_checkset( MEEG ); %notch filter around 50hz
    MEEG = pop_resample(MEEG,256);
    chans2inc = [1:1:62]; MEEG = pop_runica(MEEG,'runica','chanind',chans2inc);
    
    %==========================================================================

    pop_saveset(MEEG,'filename',['T1_' currentSubject '_SA_allbrush_0504.set'],'filepath',currentDirectory);

    disp(['Subject T' currentSubject ' epoched, merged file saved after some preprocessing'])
end
