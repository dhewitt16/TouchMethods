%This program allows viewing of grand average TF data from the brushing
%experiment. It uploads four conditions (slow palm, slow arm, fast palm,
%fast arm), applies the ERD formula to the power data, and creates
%topoplots and time frequency plots showing grand average results.

%Originally by Andrej Stancak, 20/11/2020
%Last modified by Danielle Hewitt, 11/04/23

%=====================================================================
% Set paths for EEGLAB AND FIELDTRIP
eeglab_path = '/Users/dhewitt/Analysis/eeglab2023.1/'; addpath(eeglab_path);
fieldtrip_path = '/Users/dhewitt/Analysis/fieldtrip-20240110/'; addpath(fieldtrip_path); ft_defaults;
d = char(datetime('today'));
%=====================================================================
close all;

cfg  = [];
cfg.dir  = '/Users/dhewitt/OneDrive - Nexus365/Data/Touch/UK/TouchStudy1/Preprocessed_0905/';
cfg.sub  = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15'};
%cfg.sub  = {'01'};
cfg.conds = {'sPalm', 'sArm', 'fPalm', 'fArm'};

alpha = [8 12]; beta = [16 24]; theta = [4 7]; %setting frequency bands of interest
prefrontal = {'Fp1','Fp2','AFz'}; frontal = {'F1', 'F2', 'F3', 'F4'}; central = {'C3', 'C1', 'Cz', 'C2', 'C4'}; parietal = {'P3', 'P1', 'Pz', 'P2', 'P4'}; temporal = {'T7', 'T8'}; occipital = {'O1','O2'}; %optional electrode clusters

cfg.fband = theta; %change to determine which band to run on = run for each band of interest

%=====================================================================
%%Statistics and figures
stats = 0; %for permutation and excel file export, change to 1. otherwise, change to 0
cfg.el = {'C3','C1','CP3','P5','P3','P1'}; %type in electrodes themselves or choose a cluster
cfg.epochTime = [0.5 3]; %for stats
cfg.baselineTime = [-2 -1]; %for baseline correction
cfg.plotTime = [-2 4.9];
cfg.scale = 30; %scale for figures

%=====================================================================

disp(['This analysis involves ' num2str(size(cfg.sub,2)) ' participants']);

%loads .mat containing D struct with all pts data if script has already been
%run. Also loads Gfreq and Gtime arrays
Gname = [cfg.dir 'Touch_allspectra_020524.mat'];

if exist(Gname)==2
    load(Gname);
    disp(['File ' Gname ' with ' num2str(size(D.fastArm,4)) ' has been loaded']);
    G.bfArm.time = Gtime;
    G.bfArm.freq = Gfreq;
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

    for i=1:size(cfg.sub,2) %%doing the same data twice/for each participant - check why

       % wname = ([cfg.dir 'T1_' cfg.sub{i} '/' 'T1_' cfg.sub{i} '_allbrush_3004.mat']);
        wname = ([cfg.dir 'T1_' cfg.sub{i} '_allbrush_3004.mat']);
        if exist(wname)==0
            disp(['non-existent file: ' wname]);
            return;
        end

        load(wname); %load each cfg.subject until all done .
        disp([' just loaded  ' wname]);

        msPalm = G.ERDsPalm;
        D.slowPalm(:,:,:,i)=msPalm(:,:,:);

        mBsArm = G.ERDsArm;
        D.slowArm(:,:,:,i)=mBsArm(:,:,:);

        mBfPalm = G.ERDfPalm;
        D.fastPalm(:,:,:,i)=mBfPalm(:,:,:);

        mBfArm = G.ERDfArm;
        D.fastArm(:,:,:,i)=mBfArm(:,:,:);

    end

    save(Gname,'D');
    Gfreq= G.BfArm.freq;
    save(Gname,'Gfreq','-Append');
    Gtime= G.BfArm.time;
    save(Gname,'Gtime','-Append');

    disp(['File ' Gname ' with ' num2str(size(cfg.sub,2)) ' has been saved for future analyses']);

end

disp('All files read ok');
%=========================================================================

E = readlocs('/Users/dhewitt/Analysis/TouchStudy/chanlocs.sfp');

s1=nearest(Gtime,cfg.epochTime(1)); %getting time windows
s2=nearest(Gtime,cfg.epochTime(2));

f1=nearest(Gfreq,cfg.fband(1)); %getting freq band
f2=nearest(Gfreq,cfg.fband(2));

VslowPalm = squeeze(mean(D.slowPalm(:,f1:f2,:,:),2)); %averaging over freq band
VslowPalm = squeeze(mean(VslowPalm(:,s1:s2,:),2)); %averaging over time window - now 63 els x pts

VslowArm = squeeze(mean(D.slowArm(:,f1:f2,:,:),2));
VslowArm = squeeze(mean(VslowArm(:,s1:s2,:),2));

VfastPalm = squeeze(mean(D.fastPalm(:,f1:f2,:,:),2));
VfastPalm = squeeze(mean(VfastPalm(:,s1:s2,:),2));

VfastArm = squeeze(mean(D.fastArm(:,f1:f2,:,:),2));
VfastArm = squeeze(mean(VfastArm(:,s1:s2,:),2));

%==========================================================================
%plotting the topoplots for grand averages for all conditions
close all;

figure('Name','Grand average topoplots for each condition');

subplot(2,2,1);
topoplot(mean(-VslowPalm,2),E,'style','map','maplimits',[-cfg.scale cfg.scale]); colorbar;
title('All Slow Palm');

subplot(2,2,2);
topoplot(mean(-VslowArm,2),E,'style','map','maplimits',[-cfg.scale cfg.scale]); colorbar;
title('All Slow Arm');

subplot(2,2,3);
topoplot(mean(-VfastPalm,2),E,'style','map','maplimits',[-cfg.scale cfg.scale]); colorbar;
title('All Fast Palm');

subplot(2,2,4);
topoplot(mean(-VfastArm,2),E,'style','map','maplimits',[-cfg.scale cfg.scale]); colorbar;
title('All Fast Arm');

%creating blank map showing electrode labels
figure('Name','Blank map showing electrode locations');
v=zeros(63,0);
topoplot(v,E,'electrodes','labels');

%==========================================================================
%This section extrxts the indices of all valid electrode labels in cfg.els

EL = [];
for j=1:size(cfg.el,2)
    for k=1:63
        if strcmp(cfg.el{j},E(k).labels)==1
            EL = [EL k];
        end
    end
end

%==========================================================================
%plotting the time-frequency plots for one or more electrodes in four
%conditions

%averaging over participants
M.slowPalm = mean(D.slowPalm,4);
M.slowArm = mean(D.slowArm,4);
M.fastPalm = mean(D.fastPalm,4);
M.fastArm = mean(D.fastArm,4);

s1=nearest(Gtime,cfg.plotTime(1)); %getting time windows
s2=nearest(Gtime,cfg.plotTime(2));

figure('Name', 'Time frequency plots for each condition');

CLIM = [-cfg.scale cfg.scale];
subplot(2,2,1);
imagesc(squeeze(mean(-M.slowPalm(EL,:,s1:s2),1)),CLIM); colorbar;
title('Slow Palm');

subplot(2,2,2);
imagesc(squeeze(mean(-M.slowArm(EL,:,s1:s2),1)),CLIM); colorbar;
title('Slow Arm');

subplot(2,2,3);
imagesc(squeeze(mean(-M.fastPalm(EL,:,s1:s2),1)),CLIM); colorbar;
title('Fast Palm');

subplot(2,2,4);
imagesc(squeeze(mean(-M.fastArm(EL,:,s1:s2),1)),CLIM); colorbar;
title('Fast Arm');


%% =================================================================
%export for statistics

if stats == 0
    disp('Results not saved to file');
    return
else

    % Running permutation analysis using statcond from EEGLAB
    [F df P] = statcond({VslowPalm VslowArm; VfastPalm VfastArm},'method','perm','naccu',5000); %permutation test

    % Plotting significant electrodes
    MESpeed = P{1};
    for sigel = 1:63

        if MESpeed(sigel,1) > 0.05
            MESpeed(sigel,1) = 0;
        else
            if MESpeed(sigel,1) <= 0.05
                MESpeed(sigel,1) = 100;
            end
        end

    end

    %Main effect of touch loc
    MELoc = P{2};
    for sigel = 1:63

        if MELoc(sigel,1) > 0.05
            MELoc(sigel,1) = 0;
        else
            if MELoc(sigel,1) <= 0.05
                MELoc(sigel,1) = 100;
            end
        end

    end

    %Interaction between touch location and speed
    Int = P{3};
    for sigel = 1:63

        if Int(sigel,1) > 0.05
            Int(sigel,1) = 0;
        else
            if Int(sigel,1) <= 0.05
                Int(sigel,1) = 100;
            end
        end

    end

    figure('Name', 'Permutation Results');
    subplot(2,2,1); topoplot(MESpeed,E,'style','map'); title('Main Effect Speed');
    subplot(2,2,2); topoplot(MELoc,E,'style','map'); title('Main Effect Location');
    subplot(2,2,3); topoplot(Int,E,'style','map'); title('Speed x Location');

    %=====================================================================
    %now creating the export file
    S = zeros(size(cfg.sub,2),4);
    for i=1:size(cfg.sub,2)
        m=squeeze(mean(D.slowPalm(EL,f1:f2,s1:s2,i)));
        m=squeeze(mean(m));
        m=squeeze(mean(m));
        S(i,1)=m;

        m=squeeze(mean(D.slowArm(EL,f1:f2,s1:s2,i)));
        m=squeeze(mean(m));
        m=squeeze(mean(m));
        S(i,2)=m;

        m=squeeze(mean(D.fastPalm(EL,f1:f2,s1:s2,i)));
        m=squeeze(mean(m));
        m=squeeze(mean(m));
        S(i,3)=m;

        m=squeeze(mean(D.fastArm(EL,f1:f2,s1:s2,i)));
        m=squeeze(mean(m));
        m=squeeze(mean(m));
        S(i,4)=m;

    end

    elstring = cfg.el{1};
    for k=2:size(cfg.el,2)
        elstring = [elstring '_' cfg.el{k}];
    end

    outname = [cfg.dir elstring '_' num2str(f1) '-' num2str(f2) 'Hz_' num2str(s1) '-' num2str(s2) 'ms.xlsx'];
    writematrix(S,outname);
    disp(['Results saved to ' outname])

end
%=====================================================================
