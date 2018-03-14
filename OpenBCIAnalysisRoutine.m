%% Opens OpenBCI ASCII file into MATlab
% eeg_data = csvread('OpenBCI-RAW-2017-08-12_Second.txt', 6, 1);
% eeg_data = csvread('OpenBCI-RAW-2017-08-12_EyesOpen.txt', 6, 1);
% eeg_data = csvread('OpenBCI-RAW-2017-08-12_EyesOpenMeditation.txt', 7, 1);
eeg_data = csvread('OpenBCI-RAW-2018-03-13_19-54-33.txt', 7, 1);

% Row offset is the number of rows in your txt file before the start
% of your EEG data (in the current version of the OpenMedBCI GUI, there
% are 5 commented lines before the start of the data, so the offset 
% should be 5 to make the matrix start on line 6). Column offset 
% skips the sample number column.

%% Removes anomily zero value which for some reason was interspaced every second value
eeg_data(2:2:end, :)=[]; % This interspacing is not in the orgiional ASCII but only once read into matlab.
%% Just selectign channels used in this last recording: 4 and 7
eeg_data1(:, :)=eeg_data(:,1:8);


%% Rotates matrix
eeg_data=eeg_data';
%%
save eeg_data eeg_data
%% USE EEGlab to open .m file. If this doesn't have triggers allready in it
% use the follwoing to place triggers every second:

for i = 1:314;% max = seconds
    EEG.event(1,i).type = 10;
    EEG.event(1,i).latency = i*256-256+1;
    EEG.event(1,i).urevent = 10;
end

%% Then save set file from EEGLab. 
% Note, this routine seems to work, even if you do not epoch the data
% within EEGlab. This is still worth comparing, eg FT analysis with and
% without prior epoching in EEGLab. 

% Data does not need to be epoched in eeglab, as this is done below. 
%if done in eeglab and ft, then you end up with two of every epoch.


%% This converts data strucure from EEGLAB (needs to be open) into FT
% However this is not specificaly used here.
% This is because I was not able to define trial or use events in a file
% that had not been opened with pop. Not sure why. 

%This line is just included here so as not to lose it. Not used here. 
data = eeglab2fieldtrip( EEG, 'preprocessing', 'none' );

%% Loads set file which has previously had events added. Can do with or without prior epoching in EEGLab
  load.set = '*.set*';
load.title = 'Please select files for processing';
[files.raw files.path] = uigetfile(load.set, load.title, 'MultiSelect', 'on');
cd (files.path) 
%%
    cfg = [];
    cfg.dataset = files.raw;
    cfg.trialfun = 'ft_trialfun_general'; % this is the default
    cfg.trialdef.eventtype = 'trigger';
    cfg.trialdef.eventvalue = 10;% value put in trigger events
    cfg.trialdef.prestim = 0; % in seconds. THis is long, so that I can baseline from blank period priod to trial
    cfg.trialdef.poststim = 1; % in seconds
%     cfg.channel='chan001'; % canselect only one channel otherwise, all. 
    cfg = ft_definetrial(cfg);
    
%     prepro.filt(1,1)=.5;%s5
%     prepro.filt(2,1)=45;%s5
%     cfg.lpfilter = 'yes';
    cfg.bpfilter      = 'yes'  
    cfg.bpfreq        = [2 45];
%     cfg.lpfreq = prepro.filt(2,1);
%     cfg.hpfreq = prepro.filt(1,1);
    cfg.demean     = 'yes';  
%     cfg.channel    = {'1'}; % Default: all
     
preprocessed= ft_preprocessing(cfg);
%%
cfg        = [];
cfg.method = 'fastica'; % this is the default and uses the implementation from EEGLAB
% cfg.channel = {'MEG','BIO002', 'BIO001'};
cfg.fastica.maxNumIterations=1000;
comp = ft_componentanalysis(cfg, preprocessed)
%%
comp=preprocessed;
%%
cfg          = [];
    cfg.method   = 'summary';
    % cfg.gradscale = 2e-2;
    % cfg.alim     = 2e-12; 
    RejVis1        = ft_rejectvisual(cfg,comp);
%% %% View channel by channel data for manual artifact rejection

    cfg          = [];
    cfg.method   = 'channel';
    %cfg.alim     = 1e-12; 
    cfg.gradscale = 2e-2;
    RejVis2      = ft_rejectvisual(cfg,RejVis1);
    %%
     cfg          = [];
    cfg.method   = 'trial';


   CleanData  = ft_rejectvisual(cfg,RejVis2);

%% This appears to work. Havne't taken it further at present.
    cfg=[];
%  cfg.trials         = Open.trial;
    cfg.channel='chan006';
  cfg.output         = 'pow';
  cfg.method         = 'mtmconvol';
  cfg.taper          = 'hanning';
  cfg.foi            = [5:2.5:40];
  cfg.t_ftimwin      = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi            = .25:0.025:.725;  
  cfg.keeptrials = 'yes';
%   cfg.keeptrials = 'no';
     TFRhann   = ft_freqanalysis(cfg, CleanData);
%% DOES NOT WORK YET BECASUE OF LACK OF LAYOUT INFO
% cfg = [];
% cfg.baseline     = [-0.5 -0.1]; 
% cfg.baselinetype = 'absolute'; 
% cfg.zlim         = [-3e-27 3e-27];	        
% cfg.showlabels   = 'yes';	
% % cfg.layout       = 'CTF151.lay';
% figure 
% ft_multiplotTFR(cfg, TFRhann);
%% SINGLE PLOT
cfg = [];
% cfg.baseline     = [-0.5 -0.1];
cfg.baselinetype = 'absolute';  
cfg.maskstyle    = 'saturation';	
% cfg.zlim         = [-3e-27 3e-27];	
% cfg.channel      =  'fastica001';
figure 
ft_singleplotTFR(cfg, TFRhann);
 Frames = getframe; % get the current frame and store it
% ft_singleplotTFR(cfg, FreqOpen);

%% Show movie
figure
movie(Frames,1)
%% Creatgin neightbours for cluster stats. 
%Only usefull if have full topography of channels

   cfg                 = [];
% cfg.feedback        = 'yes';
cfg.method          = 'template';
cfg.layout='elec1010.lay'
%   cfg.layout='neuromag306mag.lay';
  cfg.neighbours      = ft_prepare_neighbours(cfg);
%   cfg.neighbours      = ft_prepare_neighbours(cfg, FC_19_AVG);
  
%% Stats
% This produces no sig clusters. 
% Is this because it is only with two channels, 
% so how canthere be spatial clusters?
% Realy, I just want temporal clusters
% Can that be done?

% load FreqOpen
% load FreqOpen

% cfg = [];
fg.avgovertime = 'yes';
% cfg.channel        = {'chan001'};
cfg.channel          = {'fastica003'};
cfg.latency          = 'all';
cfg.frequency        = [11:13];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';

% cfg.correctm         = 'none';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
% cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
% cfg.parameter='powspctrm';
% prepare_neighbours determines what sensors may form clusters
% cfg_neighb.method    = 'distance';
% cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, FreqOpen);


design = zeros(1,size(FreqOpen.powspctrm,1) + size(FreqOpen.powspctrm,1));
design(1,1:size(FreqOpen.powspctrm,1)) = 1;
design(1,(size(FreqOpen.powspctrm,1)+1):(size(FreqOpen.powspctrm,1)+...
  size(FreqOpen.powspctrm,1))) = 2;

% design = zeros(1,size(MEG1.powspctrm,1) + size(MEG2.powspctrm,1));
% design(1,1:size(MEG1.powspctrm,1)) = 1;
% design(1,(size(MEG1.powspctrm,1)+1):(size(MEG1.powspctrm,1)+...
%   size(MEG2.powspctrm,1))) = 2;

% Nsub = 1; %num SUbjects
% cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
% cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
% cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
% cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number


cfg.design           = design;
cfg.ivar             = 1;
% cfg.uvar                = 2; 

[stat] = ft_freqstatistics(cfg, FreqOpen, FreqOpen);
% [stat] = ft_freqstatistics(cfg, Open(1,:), Open(1,:));
stat.prob
%% FOr plotting clusters on topography. Not used if just single channel comparisons

stat.raweffect = FreqOpen.powspctrm - FreqOpen.powspctrm;

cfg = [];
cfg.alpha  = 0.025;
cfg.parameter = 'raweffect';
cfg.zlim   = [-1e-27 1e-27];
cfg.layout = 'elec1010.lay';
ft_clusterplot(cfg, stat);

%% Creating just powspctrm and reshaping. For stats. 
% This may only work with non trials data. 
% This was to look over single values, to follow as in FT turotiral on mat
% lab ttest using avg data. 
Open=FreqOpen.powspctrm;
% size(Open);
Open=Open(:,:,6:15);% cutting blank space. 
Open=reshape(Open, [1,30]);
% plot(Open(1,:));
a=[1:30];

Open=FreqOpen.powspctrm;
size(Open);
Open=Open(:,:,6:15);
Open=reshape(Open, [1,30]);
% plot(a,Open(1,:), a, Open(1,:))

% Med1=Med.powspctrm;
% size(Med1);
% Med1=Med1(:,:,6:15);
% Med1=reshape(Med1,[2,30]);
%%
plot(a, Open(1,:), a, Open(1,:))


%% News flash, this now works with the Open and Open. 
% By work I mean, guves H=1
% This is just looking at the power values only. 
% I created Open and Open in the cell 2 below this. 
% Had to use data which did not have individual trials, eg as avg
% I guess this simulates a grandavg with many subj
% Just doing stats on powspctrm only. Is that what was done in the ft
% tutorial?
%Check.

% OpenminOpen = FreqOpen.powspctrm - FreqOpen.powspctrm;
OpenminOpen = Open(1,:) - Open(1,:);
[h,p,ci,stats] = ttest(OpenminOpen, 0, 0.05) % H0: mean = 0, alpha 0.05


%% Other FT tutorial for stats and ft stats require 
% grandavfg data or multple ptp and full topogtaphy. 
% Maybe try that with next data recording
% Will also need to re create cfg.design

%% Before ending, should try to make a ft_freqstatitics work
% I have used it above, but it does not produce H=1, when Open ALPha is
% clearly much higher than Open ALpha. 
% That would be a last step here: getting ft_freqstatitics to 
% produce meaningfull results, which in this case would be
% H=1 and some significant clusters. 
