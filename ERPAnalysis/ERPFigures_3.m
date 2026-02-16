% Script for epoching cleaned data files using EEGLAB,
% computing the CSD Transformation using FieldTrip spherical
% spline interpolation, and plotting ERPs. The data can be recomputed or
% selected from previously computed files by changing 'stats'.
% FieldTrip and EEGLAB must be in the path.
% If there are issues with plotting, open EEGLAB before running.

% Last updated: Danielle Hewitt, 2nd Jan 2024.

%% Set your analysis paths and get the data
%=====================================================================
% Set paths for EEGLAB AND FIELDTRIP - modify these for your own setup
eeglab_path = '/Users/dhewitt/Analysis/eeglab2024.2/'; addpath(eeglab_path);
fieldtrip_path = '/Users/dhewitt/Analysis/fieldtrip-20240110/'; addpath(fieldtrip_path); ft_defaults;
close all;  d = char(datetime('today'));
chanlocsfile = '/Users/dhewitt/Analysis/TouchStudy/chanlocs.sfp';
%========================================================================
%Specify subjects and main directory where data is stored
currentDirectory     = '/Users/dhewitt/Data/Touch/CE-UK-Collaboration/UK/ERPAnalysis/';
subs = {'01', '02', '04', '05', '06', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17','18','19','20','21','22','23','24','25'};
computeEpochs = 0; epoch = [-1 1]; baseline = [-200 0]; %1 = epoch the data with baseline specified - only needs to be done once unless changes required
computeCSD = 0; %to compute CSD transformation - only needs to be done once unless recomputing
mostRecentFile = 'AllERP_csd_29-Sep-2025.mat'; %to reload previously computed file, insert name here
statsExports = 1; %1 = export file with 4 columns, one for each condition
%========================================================================
% Plotting & export options
plottingTimesAll = [0 800]; %for ERP figures
plottingTimesSelected = [210 310]; %150 200 fir N140; 210 310 P300 for stats exports and topographies
selectedElectrode = {'Cz','Pz'}; % 'FC1','FCz','FC2' for N140, 'Cz', 'Pz' for P300;  multiple electrodes can be specified
%========================================================================

if computeEpochs == 1
    for iSub = 1:size(subs,2)
        currentSubject = subs{iSub};

        wname = ([currentDirectory '/ICARej/T2_' currentSubject '_icarejv5.set']);
        if exist(wname) == 0
            disp(['File ' wname ' does not exist']);
            return
        end

        %==========================================================================

        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', 'Continuous EEG Data epochs');
        eeglab redraw;
        EEG = pop_loadset(wname); EEG = eeg_checkset( EEG );
        EEG = pop_epoch( EEG, {  'B1'  }, epoch, 'newname', 'epochs', 'epochinfo', 'yes');
        EEG = pop_rmbase( EEG, baseline ,[]); EEG = eeg_checkset( EEG );
        EEG = pop_epoch( EEG, {  }, [-0.2         0.8], 'newname', 'b1_epochs_blc_epochs', 'epochinfo', 'yes'); EEG = eeg_checkset( EEG );
        pop_saveset(EEG,'filename',['T2_' currentSubject '_b1epochs.set'],'filepath',currentDirectory); %saving the preproc file for later

        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', 'Continuous EEG Data epochs');
        eeglab redraw;
        EEG = pop_loadset(wname); EEG = eeg_checkset( EEG );
        EEG = pop_epoch( EEG, {  'B2'  }, epoch, 'newname', 'epochs', 'epochinfo', 'yes');
        EEG = pop_rmbase( EEG, baseline ,[]); EEG = eeg_checkset( EEG );
        EEG = pop_epoch( EEG, {  }, [-0.2         0.8], 'newname', 'b2_epochs_blc_epochs', 'epochinfo', 'yes'); EEG = eeg_checkset( EEG );
        pop_saveset(EEG,'filename',['T2_' currentSubject '_b2epochs.set'],'filepath',currentDirectory); %saving the preproc file for later

        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', 'Continuous EEG Data epochs');
        eeglab redraw;
        EEG = pop_loadset(wname); EEG = eeg_checkset( EEG );
        EEG = pop_epoch( EEG, {  'B3'  }, epoch, 'newname', 'epochs', 'epochinfo', 'yes');
        EEG = pop_rmbase( EEG, baseline ,[]); EEG = eeg_checkset( EEG );
        EEG = pop_epoch( EEG, {  }, [-0.2         0.8], 'newname', 'b3_epochs_blc_epochs', 'epochinfo', 'yes'); EEG = eeg_checkset( EEG );
        pop_saveset(EEG,'filename',['T2_' currentSubject '_b3epochs.set'],'filepath',currentDirectory); %saving the preproc file for later

        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', 'Continuous EEG Data epochs');
        eeglab redraw;
        EEG = pop_loadset(wname); EEG = eeg_checkset( EEG );
        EEG = pop_epoch( EEG, {  'B4'  }, epoch, 'newname', 'epochs', 'epochinfo', 'yes');
        EEG = pop_rmbase( EEG, baseline ,[]); EEG = eeg_checkset( EEG );
        EEG = pop_epoch( EEG, {  }, [-0.2         0.8], 'newname', 'b4_epochs_blc_epochs', 'epochinfo', 'yes'); EEG = eeg_checkset( EEG );
        pop_saveset(EEG,'filename',['T2_' currentSubject '_b4epochs.set'],'filepath',currentDirectory); %saving the preproc file for later
    end
end

if computeCSD == 1

    CSDData = zeros(63,256,4,numel(subs)); %initialising data matrix
    cfg.layout = '/Users/dhewitt/Analysis/TouchStudy/chanlocs.sfp'; %this file must be in the working directory
    EGIlay  = ft_prepare_layout(cfg); %ft_prepare_layout loeeads or creates a 2-D layout of the channel locations.

    for iSub = 1:numel(subs)
        currentSubject = subs{iSub};
        disp(['Processing subject ' currentSubject]);

        for iCond = 1:4
            currentCond = iCond;
            wname = ([currentDirectory 'T2_' currentSubject '_b' num2str(currentCond) 'epochs.set']);

            disp(['Checking file: ' wname]);

            if exist(wname) == 0
                disp(['File ' wname ' does not exist']);
                continue;
            end

            %==========================================================================
            eEEG = pop_loadset(wname);
            %==========================================================================

            FTdata = eeglab2fieldtrip(eEEG,'preprocessing'); %This function converts EEGLAB datasets to Fieldtrip for source localization

            cfg = []; %create empty matrix
            cfg.method       = 'spline'; %The 'hjorth' method implements B. Hjort; An on-line transformation of EEG scalp potentials into orthogonal source derivation. Electroencephalography and Clinical Neurophysiology 39:526-530, 1975.
            cfg.trials       = 'all'; %'all' or a selection given as a 1xN vector (default = 'all')
            cfg.conductivity = 0.33;
            cfg.lambda = 1e-5;
            cfg.degree = 14;
            cfg.feedback = 'gui';

            [FTdata] = ft_scalpcurrentdensity(cfg, FTdata); %ft_scalpcurrentdensity computes an estimate of the SCD using the second-order derivative (the surface Laplacian) of the EEG potential distribution

            %==========================================================================

            %save([[currentDirectory 'T2_' currentSubject '_b' currentCond
            %'_csd_512.mat']], 'FTdata'); % to save single subject data

            % Update the Data matrix based on the conditions
            CSDData(:, :, iCond, iSub) = mean(cat(3, FTdata.trial{:}), 3); %saving all in matrix
            Times = (cell2mat(FTdata.time(1))).*1000; %saving times for plotting later
            %Data(:, :, :,iCond, iSub) = cat(3, FTdata.trial{:});
        end

    end

    allSubOutname = ([currentDirectory 'AllERP_csd_' d '.mat']);
    save(allSubOutname, 'CSDData');
    save(allSubOutname,'Times','-Append');
else
    load([currentDirectory mostRecentFile])
   % Times = Times.*1000;
end

%==========================================================================

%% Now plotting

E = readlocs(chanlocsfile);
t1=nearest(Times,plottingTimesAll(1)); t2=nearest(Times,plottingTimesAll(2));

v=CSDData(:,t1:t2,:,:);
v = squeeze(mean(v,2));
v = squeeze(mean(v,3));

%plotting grand average topos for each condition

figure('Name','Topo grand average for each condition, 0-0.8s');
splttl={'SlowOther','SlowSelf','FastOther','FastSelf'};
for i = 1:size(v,2)
    subplot(1,5,1);topoplot(mean(v,2),E,'maplimits',[-4 4]); colorbar; title('Mean');
    subplot(1,5,i+1); topoplot(squeeze(v(:,i)), E,'maplimits',[-4 4]); colorbar;
    title(splttl{i});
end

% %grand av over all conditions
% figure('Name','Grand av of all conditions, 0-0.8s');
% topoplot(mean(v,2),E);

%==========================================================================
%This section extracts the indices of all valid electrode labels in cfg.els

EL = [];
for j=1:size(selectedElectrode,2)
    for k=1:63
        if strcmp(selectedElectrode{j},E(k).labels)==1
            EL = [EL k];
        end
    end
end

%==========================================================================

meanV = squeeze(mean(CSDData(EL,:,:,:),1)); %mean over selected els
meanV = squeeze(mean(meanV,3)); %mean over participants
allElmeanV = squeeze(mean(mean(CSDData,4),1)); %mean over all electrodes and all participants
p1 = allElmeanV(:,1); p2 = allElmeanV(:,2); p3 = allElmeanV(:,3); p4 = allElmeanV(:,4);

%% Getting standard deviations
% stdevs = squeeze(mean(CSDData(:,:,:,:),1));
% stdp1 = std([squeeze(stdevs(:,1,:))],0,2); stdp2 = std([squeeze(stdevs(:,2,:))],0,2); stdp3 = std([squeeze(stdevs(:,3,:))],0,2); stdp4 = std([squeeze(stdevs(:,4,:))],0,2);
% stdother = stdp1-stdp3; stdself = stdp2-stdp4;

%%
t = Times;
figure('Name','ERP grand average for each condition, all electrodes');
plot(t,p1,'b'); hold on;
plot(t,p2,'r'); hold on;
plot(t,p3,'g'); hold on;
plot(t,p4,'k'); hold on;
legend('SlowOther','SlowSelf','FastOther','FastSelf'); hold on; grid on;

%%  difference waves

% figure('Name','Difference Waves, all electrodes');
% other = p1-p3;
% self = p2-p4;
% plot(t,other,'b'); hold on;
% plot(t,self,'r'); hold on;
% legend('SlowOther-FastOther','SlowSelf-FastSelf'); hold on; grid on;

%% now in a particular electrode

selp1 = meanV(:,1); selp2 = meanV(:,2); selp3 = meanV(:,3); selp4 = meanV(:,4);

figure('Name','ERP grand average for each condition, selected electrode');
plot(t,selp1,'b'); hold on;
plot(t,selp2,'r'); hold on;
plot(t,selp3,'g'); hold on;
plot(t,selp4,'k'); hold on;
legend('SlowOther','SlowSelf','FastOther','FastSelf'); hold on; grid on;

%% Now topoplots in selected time window

t1=nearest(Times,plottingTimesSelected(1)); t2=nearest(Times,plottingTimesSelected(2));

v=squeeze(mean(CSDData(:,t1:t2,:,:),2));
v = squeeze(mean(v,3));

figure('Name','Topo grand average for each condition, selected time window');
splttl={'SlowOther','SlowSelf','FastOther','FastSelf'};
for i = 1:size(v,2)
    subplot(1,5,1);topoplot(mean(v,2),E,'style', 'map','plotrad', 0.5, 'maplimits',[-8 8]); colorbar; title('Mean');
    subplot(1,5,i+1); topoplot(squeeze(v(:,i)), E,'style', 'map','plotrad', 0.5,'maplimits',[-8 8]); colorbar; 
    title(splttl{i});
end


% figure('Name','Grand av for each condition, 180-210ms'); for i = 1:size(v,2)
%     subplot(1,4,i); topoplot(squeeze(v(:,i)), E);
% end
%
% %grand av over all conditions
% figure('Name','Grand average over all conditions, 180-210ms');
% topoplot(mean(v,2),E);

if statsExports == 1
    data = zeros(23,5); %adding mean vals or adding subject ID
    %  data = zeros(23,4);
    vals=squeeze(mean(CSDData(EL,t1:t2,:,:),1));
    vals=squeeze(mean(vals,1))';

    for j=1:size(vals,2)
        %----- adding mean vals
        %data(:,1) = squeeze(mean(vals,2));
        %data(:,j+1) = vals(:,j);
        %----- not adding mean vals
        %data(:,j) = vals(:,j);
        %----- with subject ID
        data(:,1) = str2double(subs)';
        data(:,j+1) = vals(:,j);
    end

    elstring = selectedElectrode{1};
    for k=2:size(selectedElectrode,2)
        elstring = [elstring '-' selectedElectrode{k}];
    end

    outname = [currentDirectory 'AllERP_csd_' elstring '_' num2str(plottingTimesSelected(1)) '_' num2str(plottingTimesSelected(2)) '.xlsx'];
    writematrix(data,outname);
    disp(['Results saved to ' outname])

end

