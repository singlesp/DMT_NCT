%% Cluster Stats for Control Energy paper
% Chris Timmermann 2024

load('RegressorLZInterpscrubbedConvolvedAvg.mat')

%baseline correct EEG first:
bldmt = nanmean(RegDMT2(:,1:240),2);
RegDMT2 = RegDMT2 - bldmt;
RegPCB2 = RegPCB2 - nanmean(RegPCB2(:,1:240),2);

RegDMT2 = RegDMT2(:,2:839);
RegPCB2 = RegPCB2(:,2:839);

run_cluster_stats(RegDMT2,RegPCB2,1,{'EEG Signal Diversity'},[1 0 0],[0 0 0]);

%% Run stats
% First grand average but keep the trials

% Load Fieldtrip
addpath  '/Users/sps253/Documents/MATLAB/fieldtrip-20230118';
ft_defaults

load('ERPtrials.mat')

cfg = [];
cfg.keepindividual = 'yes' ;
[GADMTac] = ft_timelockgrandaverage(cfg, ERPacute{1:14});
[GAPCBac] = ft_timelockgrandaverage(cfg, ERPacute{15:28});

GADMTac.individual(:,:,839:end) = [];
GADMTac.individual(1:14,1,1:838) = RegDMT2;
GADMTac.time = 1:838;

GAPCBac.individual(:,:,839:end) = [];
GAPCBac.individual(1:14,1,1:838) = RegPCB2;
GAPCBac.time = 1:838;



% select data for analysis
xdata = GADMTac;
ydata = GAPCBac;



%% Do cluster stats

% specifies with which sensors other sensors can form clusters
cfg_neighb.method    = 'template';
cfg_neighb.layout    =   'biosemi64.lay';
cfg_neighb.template = 'biosemi64_neighb.mat';
cfg_neighb.feedback = 'no'; % Yes if you want to check the neighberhood structure
neighbours       = ft_prepare_neighbours(cfg_neighb, xdata(1));
    
cfg = [];
cfg.method = 'montecarlo';       % use the Monte Carlo Method to calculate the significance probability
cfg.statistic = 'ft_statfun_depsamplesT'; % use the independent samples T-statistic as a measure to
                               % evaluate the effect at the sample level
% cfg.latency     = [0.4 0.8];
cfg.channel = 1;

cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;         % alpha level of the sample-specific test statistic that
                               % will be used for thresholding
cfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under the
                               % permutation distribution.
% cfg.minnbchan = 2;               % minimum number of neighborhood channels that is
                               % required for a selected sample to be included
                               % in the clustering algorithm (default=0).
cfg.neighbours = neighbours;   % see below
cfg.tail = 0;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail = 0;
cfg.alpha = 0.025;               % alpha level of the permutation test
cfg.numrandomization = 5000;      % number of draws from the permutation distribution

 
Nsub = size(xdata.individual,1);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

cfg.spmversion='spm12';
ft_defaults
stat = ft_timelockstatistics(cfg,xdata,ydata);

% save stat stat

%% Find significant times

% the positive cluster
[x y] =find(stat.prob<0.025);

timesig = unique(y);

startclust = timesig(1); 
endclust = timesig(end);

% Plot Results
addpath(genpath('~/Documents/MATLAB/boundedline'))

chan=1;
   
ERPmeanDMT = squeeze(nanmean(GADMTac.individual(:,chan,:),2));
ERPmeanPCB = squeeze(nanmean(GAPCBac.individual(:,chan,:),2));

ERPsemDMT = squeeze(std(ERPmeanDMT,[],1))/sqrt(size(ERPmeanDMT,1));
ERPsemPCB = squeeze(std(ERPmeanPCB,[],1))/sqrt(size(ERPmeanPCB,1));

x = mean(ERPmeanDMT,1);
y = GADMTac.time*2-480; % convert TR to seconds
z = ERPsemDMT;

x2 = mean(ERPmeanPCB,1);
z2 =ERPsemPCB;

figure; 
[l,p] = boundedline(y, x, z, y, x2, z2,'alpha', 'transparency', 0.28);

set(l(1), 'linewidth', 2, 'color', [0.5 0 0.5]);
set(p(1), 'facecolor', [0.5 0 0.5]);
set(l(2), 'linewidth', 2, 'color', [0 0 0]);
set(p(2), 'facecolor', [0 0 0]);

allfigs = allchild(gcf);
set(gca, 'linewidth', 2,'Fontsize',20, 'Box', 'on')

set(gcf, 'color', [1 1 1],'position' ,[1000         918         787         420]);

h2 = legend('DMT', 'Placebo');

allfigs = allchild(gcf);
set(h2, 'fontsize',22, 'FontWeight', 'bold', 'Box', 'off', 'FontName', 'Arial', 'location', 'northwest', 'LineWidth', 10);
xl = xlabel('Time (seconds)', 'FontWeight', 'bold','FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',26);
yl = ylabel('LZ', 'FontWeight', 'bold', 'FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',26);

hold on
h3 =    plot([0 0], [-2 5], '--r', 'LineWidth',2.5)

yli2 = ylim;
set(get(get(h3(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % turn off legend for this mf

xp2 = [startclust*2-480 endclust*2-480 endclust*2-480 startclust*2-480];
yp2 = [yli2(1) yli2(1) 5 5];
p2=patch(xp2,yp2,'k');
set(p2,'FaceAlpha',0.1);
set(get(get(p2(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % turn off legend for this mf

export_fig(sprintf('LZ.png'),'-m2.5')

%% Do it with continuous CE metric

load('regional_continuous_CE_DMT_gsr_volnorm.mat')
load('regional_continuous_CE_PCB_gsr_volnorm.mat')

%baseline correct  first? (up to you Parker)
bldmt = nanmean(global_CE_dmt(:,1:240),2);
RegDMT2 = global_CE_dmt - bldmt;
RegPCB2 = global_CE_pcb - nanmean(global_CE_pcb(:,1:240),2);

run_cluster_stats(RegDMT2,RegPCB2,{'Global Control Energy'},[1 0 0],[0 0 0]);

%% Run stats
% First grand average but keep the rials

% Load Fieldtrip
addpath  '/Users/christophertimmermann/Documents/EEGtoolbox/fieldtrip-20240515';
ft_defaults

load('ERPtrials.mat')

cfg = [];
cfg.keepindividual = 'yes' ;
[GADMTac] = ft_timelockgrandaverage(cfg, ERPacute{1:14});
[GAPCBac] = ft_timelockgrandaverage(cfg, ERPacute{15:28});

GADMTac.individual(:,:,839:end) = [];
GADMTac.individual(1:14,1,1:838) = RegDMT2;
GADMTac.time = 1:838;

GAPCBac.individual(:,:,839:end) = [];
GAPCBac.individual(1:14,1,1:838) = RegPCB2;
GAPCBac.time = 1:838;



% select data for analysis
xdata = GADMTac;
ydata = GAPCBac;



%% Do analysis


% specifies with which sensors other sensors can form clusters
cfg_neighb.method    = 'template';
cfg_neighb.layout    =   'biosemi64.lay';
cfg_neighb.template = 'biosemi64_neighb.mat';
cfg_neighb.feedback = 'no'; % Yes if you want to check the neighberhood structure
neighbours       = ft_prepare_neighbours(cfg_neighb, xdata(1));
    
cfg = [];
cfg.method = 'montecarlo';       % use the Monte Carlo Method to calculate the significance probability
cfg.statistic = 'ft_statfun_depsamplesT'; % use the independent samples T-statistic as a measure to
                               % evaluate the effect at the sample level
% cfg.latency     = [0.4 0.8];
cfg.channel = 1;

cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;         % alpha level of the sample-specific test statistic that
                               % will be used for thresholding
cfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under the
                               % permutation distribution.
% cfg.minnbchan = 2;               % minimum number of neighborhood channels that is
                               % required for a selected sample to be included
                               % in the clustering algorithm (default=0).
cfg.neighbours = neighbours;   % see below
cfg.tail = 0;                    % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail = 0;
cfg.alpha = 0.025;               % alpha level of the permutation test
cfg.numrandomization = 5000;      % number of draws from the permutation distribution

 
Nsub = size(xdata.individual,1);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

cfg.spmversion='spm12';
ft_defaults
stat = ft_timelockstatistics(cfg,xdata,ydata);

% save stat stat

%% Find significant times

% the positive cluster
[x y] =find(stat.prob<0.025); %the 1st eas the only significant positive cluster

timesig = unique(y);

startclust = timesig(1); 
endclust = timesig(end);

%% Plot Results
addpath(genpath('~/Documents/MATLAB/boundedline'))

chan=1;

ERPmeanDMT = squeeze(nanmean(GADMTac.individual(:,chan,:),2));
ERPmeanPCB = squeeze(nanmean(GAPCBac.individual(:,chan,:),2));

ERPsemDMT = squeeze(std(ERPmeanDMT,[],1))/sqrt(size(ERPmeanDMT,1));
ERPsemPCB = squeeze(std(ERPmeanPCB,[],1))/sqrt(size(ERPmeanPCB,1));

x = mean(ERPmeanDMT,1);
y = GADMTac.time*2-480; % convert TR to seconds
z = ERPsemDMT;

x2 = mean(ERPmeanPCB,1);
z2 =ERPsemPCB;

figure; 
[l,p] = boundedline(y, x, z, y, x2, z2,'alpha', 'transparency', 0.28);

set(l(1), 'linewidth', 2, 'color', [0.5 0 0.5]);
set(p(1), 'facecolor', [0.5 0 0.5]);
set(l(2), 'linewidth', 2, 'color', [0 0 0]);
set(p(2), 'facecolor', [0 0 0]);

set(gca, 'linewidth', 2,'Fontsize',20, 'Box', 'on','ylim',[-75 150])
set(gcf, 'color', [1 1 1],'position' ,[1000         918         787         420]);
h2 = legend('DMT', 'Placebo');
set(h2, 'fontsize',22, 'FontWeight', 'bold', 'Box', 'off', 'FontName', 'Arial', 'location', 'northwest', 'LineWidth', 10);
xl = xlabel('Time (seconds)', 'FontWeight', 'bold','FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',26);
yl = ylabel('Global CE', 'FontWeight', 'bold', 'FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',26);

hold on
h3 =    plot([0 0], [-75 150], '--r', 'LineWidth',2.5);

yli2 = ylim;
set(get(get(h3(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % turn off legend for this mf

xp2 = [startclust*2-480 endclust*2-480 endclust*2-480 startclust*2-480];
yp2 = [yli2(1) yli2(1) 150 150];
p2=patch(xp2,yp2,'k');
set(p2,'FaceAlpha',0.1);
set(get(get(p2(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % turn off legend for this mf

export_fig(sprintf('Global_CE.png'),'-m2.5')