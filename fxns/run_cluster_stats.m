function [num_clus,timesig] = run_cluster_stats(DMT_TS,PCB_TS,tail,ylab,color1,color2,ExportFigYN,figname,NewFigYN,Title,fontsizes)

% tail is -1 for negative clusters or +1 for positive clusters

if nargin <4
    ylab = 'Y';
end

if nargin <6
    color1=[0 0 0];
    color2=[1 0 0];
end

if nargin<7
    ExportFigYN=0;
end

if nargin<8 
    figname = 'fig';
end

if nargin<9
    NewFigYN=1;
end

if nargin<10
    fontsizes = [20 22 26]; %tics, legend, axislabels
end

%% baseline correct? y/n

%baseline correct first:
% DMT_TS = DMT_TS - nanmean(DMT_TS(:,1:240),2);
% PCB_TS = PCB_TS - nanmean(PCB_TS(:,1:240),2);

% this needs to be done prior to fxn for eeg
% DMT_TS = DMT_TS(:,2:839);
% PCB_TS = PCB_TS(:,2:839);

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
GADMTac.individual(1:14,1,1:838) = DMT_TS;
GADMTac.time = 1:838;

GAPCBac.individual(:,:,839:end) = [];
GAPCBac.individual(1:14,1,1:838) = PCB_TS;
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

timesig={};
if tail>0
    sig_clusters = find([stat.posclusters.prob]<0.025);
    num_clus = length(sig_clusters);
    for i=1:num_clus
        timesig{i} = find(stat.posclusterslabelmat==sig_clusters(i));
    end
elseif tail<0
    sig_clusters = find([stat.negclusters.prob]<0.025);
    num_clus = length(sig_clusters);
    for i=1:num_clus
        timesig{i} = find(stat.negclusterslabelmat==sig_clusters(i));
    end
end

% Plot Results
addpath(genpath('~/Documents/MATLAB/boundedline'))

chan=1;

ERPmeanDMT = squeeze(nanmean(GADMTac.individual(:,chan,:),2));
ERPmeanPCB = squeeze(nanmean(GAPCBac.individual(:,chan,:),2));

ERPsemDMT = squeeze(std(ERPmeanDMT,[],1))/sqrt(size(ERPmeanDMT,1));
ERPsemPCB = squeeze(std(ERPmeanPCB,[],1))/sqrt(size(ERPmeanPCB,1));

x = mean(ERPmeanDMT,1);
% y = GADMTac.time*2-480; % convert TR to seconds
y = GADMTac.time;
z = ERPsemDMT;

x2 = mean(ERPmeanPCB,1);
z2 =ERPsemPCB;

if NewFigYN==1
    figure;
end
[l,p] = boundedline(y, x, z, y, x2, z2,'alpha', 'transparency', 0.28);

set(l(1), 'linewidth', 2, 'color', color1);
set(p(1), 'facecolor', color1);
set(l(2), 'linewidth', 2, 'color', color2);
set(p(2), 'facecolor', color2);

allfigs = allchild(gcf);
set(gca, 'linewidth', 2,'Fontsize',fontsizes(1), 'Box', 'on')

set(gcf, 'color', [1 1 1],'position' ,[1000         918         787         420]);

h2 = legend('DMT', 'Placebo');
tics = linspace(0,28,8);
tics = tics*30;
tics(end)=838;


allfigs = allchild(gcf);
set(h2, 'fontsize',fontsizes(2), 'FontWeight', 'bold', 'Box', 'off', 'FontName', 'Arial', 'location', 'northwest', 'LineWidth', 10);
% xl = xlabel('Time (seconds)', 'FontWeight', 'bold','FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',26);
xticks(tics)
xticklabels([{'-8'},{'-4'},{'0'},{'4'},{'8'},{'12'},{'16'},{'20'}]);
xlabel('Minutes', 'FontWeight', 'bold','FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',fontsizes(3))
xlim([0 838])
yl = ylabel(ylab, 'FontWeight', 'bold', 'FontName', 'Arial', 'FontWeight', 'normal', 'fontsize',fontsizes(3));
ax = gca;
ax.Box = 'off';

if exist('Title','var')==1
    title(Title);
end

hold on
% h3 =    plot([0 0], [-2 5], '--r', 'LineWidth',2.5)

yli2 = ylim;
% set(get(get(h3(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % turn off legend for this mf

% xp2 = [startclust*2-480 endclust*2-480 endclust*2-480 startclust*2-480];
if ~isempty(timesig)
    for i=1:num_clus
        clusi = timesig{i};
        startclust = clusi(1);
        endclust = clusi(end);
        xp2 = [startclust endclust endclust startclust];
        yp2 = [yli2(1) yli2(1) yli2(2) yli2(2)];
        p2=patch(xp2,yp2,'k');
        set(p2,'FaceAlpha',0.1);
        set(get(get(p2(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % turn off legend for this mf
    end
end

if ExportFigYN==1
    addpath(genpath('~/Documents/MATLAB/export_fig'));
    export_fig(sprintf(['~/Downloads/',figname]),'-eps','-png','-transparent')
end

end