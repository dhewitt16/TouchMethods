%This program will compute the power spectral densities of EEg in 1-sec
%intervals over the period of brushing (0-5.5 s) and rest(-2.5 to 0 s). The spectra are
%computed in sliding manner with 0.9 s overlap. The input data are from
%BrainProducts system, however, preprocessed using EEGlab (filtering 1-70
%Hz, common average reference, downsampled, ICA-pruned, visual rejection of bad epochs).

close all;

%PLEASE RESTORE YOUR PATHS
%=====================================================================
% Set paths for EEGLAB AND FIELDTRIP
eeglab_path = '/Users/dhewitt/Analysis/eeglab2023.1/'; addpath(eeglab_path);
fieldtrip_path = '/Users/dhewitt/Analysis/fieldtrip-20240110/'; addpath(fieldtrip_path); ft_defaults;
%=====================================================================

%Input: a set file which needs to be changed in cfg.name

%sub  = {'01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15'};
sub  = {'06'};
dir  = '/Users/dhewitt/OneDrive - Nexus365/Data/Touch/UK/TouchStudy1/Preprocessed_0905/';
%dir = [pwd '\']; %please delete as appropriate
conds = {'sPalm', 'sArm', 'fPalm', 'fArm'};

%=====================================================================
epochStart = -2.5;
epochSteps = 0.1;
epochEnd   = 5.5;
baselineStart = -2.0; %for ERD
baselineEnd   = -1.0;
%=====================================================================

%optional plotting of data to check ERD data
fplot          = []; %figure plot structure
fplot.draw     = 1; %if 1, figures will be plotted, else only saves spectra
fplot.els      = {'C3' 'C4'};
fplot.what     = 1; %1 = ERD; 2 = power; 3 = amplitude (square root of power) - only option 1 works at the moment!
fplot.bands    = [8 13; 16 32];
fplot.erdscale = [-60 60];
%fplot.power    = [0.2 1.2]; #not used
fplot.timepnt  = [1.0 2.0 4.0]; %plotting ERD maps at 3 selected time points
fplot.clusters = {'C3' 'C5' 'C1' 'CP3' 'FC3'; 'C4' 'C2' 'C6' 'CP4' 'FC4'};

%=====================================================================

for iSub=1:size(sub,2)
    currentSubject = sub{iSub};

    wname = [dir 'T1_' currentSubject '_allbrush_ica_clean_300424-2.set'];  %pls change the backslash
    outname = [dir 'T1_' currentSubject '_allbrush_3004_2.mat'];
   % wname = [dir 'T1_' currentSubject '/' 'T1_' currentSubject '_allbrush_ica_clean.set'];  %pls change the backslash
    %outname = [dir 'T1_' currentSubject '/' 'T1_' currentSubject '_allbrush_3004.mat'];

    if exist(wname) == 0
        disp(['File : ' wname  ' does not exist']);
        return;
    end
    %=====================================================================
    %Loading data and events into FT

    eEEG = pop_loadset(wname);
    FTdata = eeglab2fieldtrip(eEEG,'preprocessing','none');

    %=====================================================================
    %Running 4 time-frequency analyses on epoched data - first selecting epochs
    %which are relevant for a block

    for i=1:4

        T = [];
        for k=1:size(eEEG.event,2)
            if strcmp(['B' num2str(i)],eEEG.event(k).type)==1
                T = [T k];
            end
        end

        cfg         = [];
        cfg.trials  =  T;

        segdata = ft_redefinetrial(cfg,FTdata);

        disp('Computing spectra...');

        %this bit gets the power within each frequency band ...

        cfg = [];
        cfg.output         = 'pow';
        cfg.channel        = 'EEG';
        cfg.method         = 'mtmconvol';
        cfg.taper          = 'dpss'; %multi-tapers - could change to hann
        cfg.keeptrials     = 'yes';
        cfg.pad            = 'nextpow2';
        cfg.foi            = [1:1:70]; %frequencies of interest changed from 1 to 100 - change back if needed
        cfg.t_ftimwin      = ones(length(cfg.foi),1).*1;
        cfg.tapsmofrq      = ones(length(cfg.foi),1).*2;
        cfg.toi            = [epochStart:epochSteps:epochEnd];
        TFRhann            = ft_freqanalysis(cfg, segdata);
        %TFRhann.powspctrm(~isfinite(TFRhann.powspctrm))=0;

        eval(['B' conds{i} '=TFRhann;']);
        TFR=TFRhann;

        %this is the ERD bit ... I guess, this could also be deleted, only the %averaging part of this code effective here

        cfg                = [];
        cfg.baseline 	   = [baselineStart baselineEnd];
        cfg.baselinetype   = 'no';  %no standardisation at this stage - we want to see absolute power values first, ERD can be computed anytime later
        TFR.dimord         = 'chan_freq_time';
        TFR.powspctrm      = squeeze(mean(TFR.powspctrm,1));
        TFR.cumtapcnt      = squeeze(TFR.cumtapcnt(1,:));
        %TFRrbl             = ft_freqbaseline(cfg,TFR);
        %TFRrbl.powspctrm   = TFRrbl.powspctrm*100; %no: I would store the absolute values of power (averaged across trials) and compute ERD separately

        eval(['B' conds{i} 'rbl=TFR;']);

    end

    G = [];
    G.BsPalm  = BsPalmrbl;  %slow palm
    G.BsArm   = BsArmrbl;   %slow arm
    G.BfPalm  = BfPalmrbl;  %fast palm
    G.BfArm   = BfArmrbl;   %fast arm

    %computing ERD now

    %getting the time bins correcsponding to the resting interval
    bt1 = nearest(TFR.time,baselineStart); bt2 = nearest(TFR.time,baselineEnd);
    for i=1:4

        eval (['D = G.B' conds{i} ';']);
        RestVals = mean(D.powspctrm(:,:,bt1:bt2),3);
        erd = zeros(size(D.powspctrm,1),size(D.powspctrm,2),size(D.powspctrm,3));
        for els=1:63
            for freqs=1:size(D.powspctrm,2)
                for tms = 1:size(D.powspctrm,3)
                    erd(els,freqs,tms) = 100*(RestVals(els,freqs)-D.powspctrm(els,freqs,tms))./RestVals(els,freqs); %this is the Pfu's formula 100*(R-A)/R
                end
            end
        end

        eval(['G.ERD' conds{i} '=erd ;']);

    end

    save(outname,'G');
    disp(['All time frequency analyses were saved ok: ' outname]);

end

%%Plotting the optional figures of ERD/power
if fplot.draw==1

    %getting the electrode indices
    elecs = [];
    for i=1:2
        for k=1:63
            if strcmp(D.label{k},fplot.els(1,i))==1
                elecs=[elecs k];
            end
        end
    end

    %frequency bands indices
    inds = [];
    inds(1,1) = nearest(D.freq,fplot.bands(1,1));
    inds(1,2) = nearest(D.freq,fplot.bands(1,2));
    inds(2,1) = nearest(D.freq,fplot.bands(2,1));
    inds(2,2) = nearest(D.freq,fplot.bands(2,2));

    %Figure 1 plots ERD curves in 2 electrodes and 2 bands
    figure('Name','ERD curves');
    c=0;
    for band = 1:2
        for i=1:4 %conditions

            c=c+1;
            subplot(2,4,c);
            eval(['erd = G.ERD' conds{i} ';']);
            e1 = squeeze(squeeze(mean(erd(elecs(1,1),inds(band,1):inds(band,2),:),2)));
            e2 = squeeze(squeeze(mean(erd(elecs(1,2),inds(band,1):inds(band,2),:),2)));
            plot(D.time,-e1,'r'); %while ERD is a power decrease, it is actualy a positive value; therefore, for the purpose of plotting, the ERD curve is actually plotted withinverse sign (not logigal, huh?)
            hold on;
            plot(D.time,-e2,'b');
            legend(fplot.els{1},fplot.els{2},'location','best');
            axis([-2 4.9 fplot.erdscale(1,1) fplot.erdscale(1,2)]);
            grid on;
            title([conds{i} '  ' num2str(inds(band,1)) '-' num2str(inds(band,2)) ' Hz']);

        end
    end

    %Figures 2 and 3 plot the 2D maps

    E = readlocs('chanlocs.sfp'); %please check the presence of the electrode coordinate file when calling rlocs64

    for band = 1:2
        figure('Name',['Maps1 in ' num2str(inds(band,1)) '-' num2str(inds(band,2))]); %each band on a new figure

        c=0;
        for k=1:4
            eval(['erd = G.ERD' conds{k} ';']);
            for i=1:3
                c=c+1;
                subplot(2,6,c);
                t1=nearest(D.time,fplot.timepnt(i));
                v = squeeze(squeeze(mean(erd(:,inds(band,1):inds(band,2),t1),2)));
                topoplot(-v,E,'style','map','maplimits',fplot.erdscale);  %the maps are the relative power maps, not ERD maps P[%] = -ERD[%]
                title([conds{k} ' ' num2str(fplot.timepnt(i)) ' s']);
            end
        end

    end


    %finally, time-frequency plots from clusters of electrodes around C3 and C4
    %in each of 4 conditions, defined in fplot.clusters

    %finding the electrode indices for the clusters of electrodes
    clind = zeros(5,2);
    for i=1:5
        for j=1:2
            for k=1:63
                if strcmp(D.label{k},char(fplot.clusters{j,i}))==1
                    clind(i,j)=k;
                end
            end
        end
    end


    figure('Name','TFR plots');
    c=0;
    for k=1:4 %4 conditions
        eval(['erd = G.ERD' conds{k} ';']);
        for cl=1:2
            c=c+1;
            tfr=squeeze(mean(erd(squeeze(clind(cl,:)),:,:),1));
            subplot(4,2,c);
            imagesc(tfr,fplot.erdscale);
            ax=gca;
            ax.XLim = [11 70];
            ax.YLim = [1 70];
            set(gca,'XLim',ax.XLim,'YLim',ax.YLim);
            title([conds{k} ' ' fplot.clusters{cl,1}]);
            colorbar;
        end

    end

end %end of fplot option
