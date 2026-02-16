
%% Script to export all Touch ERD data as a table in wide format.

% The input is data from all participants, which can be recomputed here or
% loaded from file.

% Output is a table with 9 columns:
% Participants x FreqBand (first value 4, 8 or 16) x FreqBandName x
% Electrode (63) x ElectrodeName x SlowOther x SlowSelf x FastOther x
% FastSelf

% Danielle Hewitt, 12th Jan 2024

%=====================================================================
% Set paths for EEGLAB AND FIELDTRIP
eeglab_path = '/Users/dhewitt/Analysis/eeglab2023.1/'; addpath(eeglab_path);
fieldtrip_path = '/Users/dhewitt/Analysis/fieldtrip-20240110/'; addpath(fieldtrip_path); ft_defaults;
d = char(datetime('today'));
%=====================================================================

cfg  = [];
cfg.sub = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15'};
cfg.dir  = '/Users/dhewitt/OneDrive - Nexus365/Data/Touch/UK/TouchStudy1/Preprocessed_0905/';
cfg.epochTime = [0.5 3]; %for stats
cfg.fband    = [4 7; 8 13; 16 24];
%=====================================================================

disp(['This analysis involves ' num2str(size(cfg.sub,2)) ' participants']);

%loads .mat containing D struct with all pts data if script has already been
%run. Also loads Gfreq and Gtime arrays
Gname = [cfg.dir 'Touch_allspectra_020524.mat'];

if exist(Gname)==2
    load(Gname);
    disp(['File ' Gname ' with ' num2str(size(D.fastArm,4)) ' has been loaded']);
    % G.bfArm.time = Gtime;
    % G.bfArm.freq = Gfreq;
    if size(D.fastArm,4)< size(cfg.sub,2)
        disp('Warning !!! there seem to be larger number of participants than in the uploaded file');
        return
    end

else %if the .mat file does not exist, creates a new D struct for the averaged ERD data to be placed in

    D.slowPalm = zeros(63,70,81,size(cfg.sub,2));
    D.slowArm = zeros(63,70,81,size(cfg.sub,2));
    D.fastPalm = zeros(63,70,81,size(cfg.sub,2));
    D.fastArm = zeros(63,70,81,size(cfg.sub,2));

    disp('wait, loading spectra');

    for i=1:size(cfg.sub,2)

        wname = ([cfg.dir 'T1_' cfg.sub{i} '/T1_' cfg.sub{i} '_allbrush_0405.mat']);
        if exist(wname)==0
            disp(['non-existent file: ' wname]);
            return;
        end

        load(wname); %load each cfg.subject until all done .
        disp([' just loaded  ' wname]);

        msPalm = G.ERDsPalm; D.slowPalm(:,:,:,i)=msPalm(:,:,:);
        mBsArm = G.ERDsArm; D.slowArm(:,:,:,i)=mBsArm(:,:,:);
        mBfPalm = G.ERDfPalm; D.fastPalm(:,:,:,i)=mBfPalm(:,:,:);
        mBfArm = G.ERDfArm; D.fastArm(:,:,:,i)=mBfArm(:,:,:);

    end

    save(Gname,'D');
    Gfreq= G.BfArm.freq; save(Gname,'Gfreq','-Append');
    Gtime= G.BfArm.time; save(Gname,'Gtime','-Append');
    disp(['File ' Gname ' with ' num2str(size(cfg.sub,2)) ' has been saved for future analyses']);

end

disp('All files read ok');
%=========================================================================

bigexport = zeros([numel(cfg.sub)*63*size(cfg.fband,1),7]);
c = 0;

for freqi = 1:size(cfg.fband,1)
    freqBand = cfg.fband(freqi,:);

    s1=nearest(Gtime,cfg.epochTime(1)); s2=nearest(Gtime,cfg.epochTime(2));
    f1=nearest(Gfreq,freqBand(1)); f2=nearest(Gfreq,freqBand(2));

    VslowPalm = squeeze(mean(D.slowPalm(:,f1:f2,:,:),2)); %averaging over freq band
    VslowPalm = squeeze(mean(VslowPalm(:,s1:s2,:),2)); %averaging over time window - now 63 els x pts
    VslowPalm=VslowPalm';
    VslowPalm=reshape(VslowPalm,[],1);

    VslowArm = squeeze(mean(D.slowArm(:,f1:f2,:,:),2));
    VslowArm = squeeze(mean(VslowArm(:,s1:s2,:),2));
    VslowArm=VslowArm';
    VslowArm=reshape(VslowArm,[],1);

    VfastPalm = squeeze(mean(D.fastPalm(:,f1:f2,:,:),2));
    VfastPalm = squeeze(mean(VfastPalm(:,s1:s2,:),2));
    VfastPalm=VfastPalm';
    VfastPalm=reshape(VfastPalm,[],1);

    VfastArm = squeeze(mean(D.fastArm(:,f1:f2,:,:),2));
    VfastArm = squeeze(mean(VfastArm(:,s1:s2,:),2));
    VfastArm=VfastArm';
    VfastArm=reshape(VfastArm,[],1);

    allconds = [VslowPalm, VslowArm, VfastPalm, VfastArm];

    %=========================================================================

    participantIDs = str2double(cfg.sub);
    participantIDs=repmat(participantIDs,1,63)';
    elIDs = [1:63];
    elIDs=repelem(elIDs,numel(cfg.sub))';

    allconds(:,5)=participantIDs;
    allconds(:,6)=freqBand(1);
    allconds(:,7)=elIDs;

    bigexport(c + (1:size(allconds,1)), :) = allconds;
    c = c + size(allconds,1);  % Update the index

end

%=========================================================================

dataTable = array2table(bigexport, 'VariableNames', {'SlowPalm', 'SlowArm', 'FastPalm', 'FastArm','ParticipantID','FreqBand','Electrode'});

elID = dataTable.Electrode;
electrodeNames={'Fp1';'F3';'F7';'FT9';'FC5';'FC1';'C3';'T7';'TP9';'CP5';'CP1';'Pz';'P3';'P7';'O1';'Oz';'O2';'P4';'P8';'TP10';'CP6';'CP2';'Cz';'C4';'T8';'FT10';'FC6';'FC2';'F4';'F8';'Fp2';'AF7';'AF3';'AFz';'F1';'F5';'FT7';'FC3';'C1';'C5';'TP7';'CP3';'P1';'P5';'PO7';'PO3';'POz';'PO4';'PO8';'P6';'P2';'CPz';'CP4';'TP8';'C6';'C2';'FC4';'FT8';'F6';'AF8';'AF4';'F2';'FCz'};
% electrodeName = electrodeNames(elID); % Map 'elID' to 'electrodeName'
% dataTable.ElectrodeName = electrodeName; %% this isn't working anymore
% for some reason

electrodeIndices = 1:numel(electrodeNames);
electrodeMap = containers.Map(electrodeIndices, electrodeNames);
getElectrodeName = @(index) electrodeMap(index);
dataTable.ElectrodeName = arrayfun(getElectrodeName, dataTable.Electrode, 'UniformOutput', false);

FreqBand = dataTable.FreqBand;
FreqBandName = cell(size(FreqBand)); % Map 'freqBand' to 'FreqBandName'
for i = 1:length(FreqBand)
    if FreqBand(i) == 4
        FreqBandName{i} = 'theta';
    elseif FreqBand(i) == 8
        FreqBandName{i} = 'alpha';
    elseif FreqBand(i) == 16
        FreqBandName{i} = 'beta';
    else
        FreqBandName{i} = 'other';
    end
end
dataTable.FreqBandName = FreqBandName;

dataTable = dataTable(:, {'ParticipantID', 'FreqBand', 'FreqBandName', 'Electrode', 'ElectrodeName', 'SlowPalm', 'SlowArm', 'FastPalm', 'FastArm'});
dataTable.Properties.VariableNames = strrep(dataTable.Properties.VariableNames, '''', '');
dataTable.ElectrodeName = strrep(dataTable.ElectrodeName, '''', '');

%=========================================================================
outname = [cfg.dir 'AllERD_UK1_' d '.xlsx'];
%outname = [cfg.dir 'allERD-allchans.csv'];
writetable(dataTable,outname);
disp(['Results saved to ' outname])

%=========================================================================



