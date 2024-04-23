%% dominance analysis - need to run min8_regional_CE_prepost_diff.m and time_resolved_CE_metrics.m first
clear all; close all;

basedir = '~/Documents/GIT/DMT_NCT/';

note='_gsr_volnorm'; %track different proc streams

diffE_dom = readtable([basedir,'results/dominance/diffE_continuous',note,'_DominanceStats.csv']);
int_dom = readtable([basedir,'results/dominance/intcorr',note,'_DominanceStats.csv']);
lz_dom = readtable([basedir,'results/dominance/LZcorr',note,'_DominanceStats.csv']);

receptors = [{'5-HT2a'},{'5-HT1a'},{'5-HT1b'},{'5-HT4'},{'5-HTT'}];
varnames=[{'HT2a'},{'HT1a'},{'HT1b'},{'HT4'},{'HTT'}];

for i = 1:length(receptors)
   
    diffE_relImp(i) = diffE_dom.PercentageRelativeImportance(strcmp(diffE_dom.Var1,varnames(i)));
    int_relImp(i) = int_dom.PercentageRelativeImportance(strcmp(int_dom.Var1,varnames(i)));
    lz_relimp(i) = lz_dom.PercentageRelativeImportance(strcmp(lz_dom.Var1,varnames(i)));
end

%% FIGURE 4

P = [diffE_relImp; int_relImp; lz_relimp];

% Spider plot
figure;
spider_plot(P,...
    'AxesLimits', [0, 0, 0, 0, 0; 50, 50, 50, 50, 50],...
    'AxesLabels',receptors);

% Legend properties
legend('control energy', 'corr(CE, intensity)', 'corr(CE, LZ)', 'Location', 'southoutside');

%% PLACEBO Dom

clear all;

basedir = '~/Documents/GIT/DMT_NCT/';

note='_gsr_volnorm'; %track different proc streams

diffE_dom = readtable([basedir,'results/dominance/diffE_pcb_continuous',note,'_DominanceStats.csv']);
int_dom = readtable([basedir,'results/dominance/intcorr_pcb',note,'_DominanceStats.csv']);
lz_dom = readtable([basedir,'results/dominance/LZcorr_pcb',note,'_DominanceStats.csv']);

receptors = [{'5-HT2a'},{'5-HT1a'},{'5-HT1b'},{'5-HT4'},{'5-HTT'}];
varnames=[{'HT2a'},{'HT1a'},{'HT1b'},{'HT4'},{'HTT'}];

for i = 1:length(receptors)
   
    diffE_relImp(i) = diffE_dom.PercentageRelativeImportance(strcmp(diffE_dom.Var1,varnames(i)));
    int_relImp(i) = int_dom.PercentageRelativeImportance(strcmp(int_dom.Var1,varnames(i)));
    lz_relimp(i) = lz_dom.PercentageRelativeImportance(strcmp(lz_dom.Var1,varnames(i)));
end

%% FIGURE 4

P = [diffE_relImp; int_relImp; lz_relimp];

% Spider plot
figure;
spider_plot(P,...
    'AxesLimits', [0, 0, 0, 0, 0; 50, 50, 50, 50, 50],...
    'AxesLabels',receptors);

% Legend properties
legend('control energy', 'corr(CE, intensity)', 'corr(CE, LZ)', 'Location', 'southoutside');