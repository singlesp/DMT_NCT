%% plot global and network level CE - must run time_resolved_CE_metrics.m first

clear all; close all;

nsub=14;

basedir = '~/Documents/GIT/DMT_NCT/';

note='_gsr_volnorm'; %track different proc streams

load([basedir,'results/regional_continuous_CE_DMT',note,'.mat'])
load([basedir,'results/regional_continuous_CE_PCB',note,'.mat'])

load([basedir,'data/sch116_to_yeo.csv'])
network7labels=sch116_to_yeo;
networks = [{'VIS'},{'SOM'},{'DAT'},{'VAT'},{'LIM'},{'FPN'},{'DMN'},{'SUB'}];

%% get network level CE's

vis_CE = sum(regional_CE_dmt(:,:,network7labels==1),3);
som_CE = sum(regional_CE_dmt(:,:,network7labels==2),3);
dat_CE = sum(regional_CE_dmt(:,:,network7labels==3),3);
vat_CE = sum(regional_CE_dmt(:,:,network7labels==4),3);
lim_CE = sum(regional_CE_dmt(:,:,network7labels==5),3);
fpn_CE = sum(regional_CE_dmt(:,:,network7labels==6),3);
dmn_CE = sum(regional_CE_dmt(:,:,network7labels==7),3);
sub_CE = sum(regional_CE_dmt(:,:,network7labels==8),3);

vis_CE_pcb = sum(regional_CE_pcb(:,:,network7labels==1),3);
som_CE_pcb = sum(regional_CE_pcb(:,:,network7labels==2),3);
dat_CE_pcb = sum(regional_CE_pcb(:,:,network7labels==3),3);
vat_CE_pcb = sum(regional_CE_pcb(:,:,network7labels==4),3);
lim_CE_pcb = sum(regional_CE_pcb(:,:,network7labels==5),3);
fpn_CE_pcb = sum(regional_CE_pcb(:,:,network7labels==6),3);
dmn_CE_pcb = sum(regional_CE_pcb(:,:,network7labels==7),3);
sub_CE_pcb = sum(regional_CE_pcb(:,:,network7labels==8),3);

figure;
subplot(2,1,1)
hold on
    plot(mean(vis_CE)','LineWidth',2)
    plot(mean(som_CE)','LineWidth',2)
    plot(mean(dat_CE)','LineWidth',2)
    plot(mean(vat_CE)','LineWidth',2)
    plot(mean(lim_CE)','LineWidth',2)
    plot(mean(fpn_CE)','LineWidth',2)
    plot(mean(dmn_CE)','LineWidth',2)
    plot(mean(sub_CE)','LineWidth',2)
legend(networks);



subplot(2,1,2)
hold on
    plot(mean(vis_CE_pcb)','LineWidth',2)
    plot(mean(som_CE_pcb)','LineWidth',2)
    plot(mean(dat_CE_pcb)','LineWidth',2)
    plot(mean(vat_CE_pcb)','LineWidth',2)
    plot(mean(lim_CE_pcb)','LineWidth',2)
    plot(mean(fpn_CE_pcb)','LineWidth',2)
    plot(mean(dmn_CE_pcb)','LineWidth',2)
    plot(mean(sub_CE_pcb)','LineWidth',2)
legend(networks);

%% SI FIGURE 1 

tics = linspace(0,28,15);
tics = tics*30;
tics(end)=838;


figure;
subplot(5,2,1.5)
[h,p]=ttest(global_CE_dmt,global_CE_pcb);
% idx = double(find(p<0.05));
post_inj = p(239:end);
pfdr = mafdr(post_inj,'BH',1);
pfdr = [ones(1,238) pfdr];
idx = double(find(pfdr<0.05));
hold on
    plot(mean(global_CE_dmt)','black','LineWidth',2)
    plot(mean(global_CE_pcb)','LineWidth',2)
text(idx',repelem(.9*max(mean(global_CE_pcb)),length(idx)),'*')
title('Global')
legend('DMT','PCB')
xticks(tics)
xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
xlabel('Minutes')


subplot(5,2,3)
[h,p]=ttest(vis_CE,vis_CE_pcb);
% idx = double(find(p<0.05));
post_inj = p(239:end);
pfdr = mafdr(post_inj,'BH',1);
pfdr = [ones(1,238) pfdr];
idx = double(find(pfdr<0.05));
hold on
    plot(mean(vis_CE)','black','LineWidth',2)
    plot(mean(vis_CE_pcb)','LineWidth',2)
text(idx',repelem(.9*max(mean(vis_CE_pcb)),length(idx)),'*')
title(networks{1})
% legend('DMT','PCB')
xticks(tics)
xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
xlabel('Minutes')


subplot(5,2,4)
[h,p]=ttest(som_CE,som_CE_pcb);
% idx = double(find(p<0.05));
post_inj = p(239:end);
pfdr = mafdr(post_inj,'BH',1);
pfdr = [ones(1,238) pfdr];
idx = double(find(pfdr<0.05));
hold on
    plot(mean(som_CE)','black','LineWidth',2)
    plot(mean(som_CE_pcb)','LineWidth',2)
text(idx',repelem(.9*max(mean(som_CE_pcb)),length(idx)),'*')
title(networks{2})
% legend('DMT','PCB')
xticks(tics)
xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
xlabel('Minutes')


subplot(5,2,5)
[h,p]=ttest(dat_CE,dat_CE_pcb);
% idx = double(find(p<0.05));
post_inj = p(239:end);
pfdr = mafdr(post_inj,'BH',1);
pfdr = [ones(1,238) pfdr];
idx = double(find(pfdr<0.05));
hold on
    plot(mean(dat_CE)','black','LineWidth',2)
    plot(mean(dat_CE_pcb)','LineWidth',2)
text(idx',repelem(.9*max(mean(dat_CE_pcb)),length(idx)),'*')
title(networks{3})
% legend('DMT','PCB')
xticks(tics)
xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
xlabel('Minutes')


subplot(5,2,6)
[h,p]=ttest(vat_CE,vat_CE_pcb);
% idx = double(find(p<0.05));
post_inj = p(239:end);
pfdr = mafdr(post_inj,'BH',1);
pfdr = [ones(1,238) pfdr];
idx = double(find(pfdr<0.05));
hold on
    plot(mean(vat_CE)','black','LineWidth',2)
    plot(mean(vat_CE_pcb)','LineWidth',2)
text(idx',repelem(.9*max(mean(vat_CE_pcb)),length(idx)),'*')
title(networks{4})
% legend('DMT','PCB')
xticks(tics)
xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
xlabel('Minutes')


subplot(5,2,7)
[h,p]=ttest(lim_CE,lim_CE_pcb);
% idx = double(find(p<0.05));
post_inj = p(239:end);
pfdr = mafdr(post_inj,'BH',1);
pfdr = [ones(1,238) pfdr];
idx = double(find(pfdr<0.05));
hold on
    plot(mean(lim_CE)','black','LineWidth',2)
    plot(mean(lim_CE_pcb)','LineWidth',2)
text(idx',repelem(.9*max(mean(lim_CE_pcb)),length(idx)),'*')
title(networks{5})
% legend('DMT','PCB')
xticks(tics)
xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
xlabel('Minutes')


subplot(5,2,8)
[h,p]=ttest(fpn_CE,fpn_CE_pcb);
% idx = double(find(p<0.05));
post_inj = p(239:end);
pfdr = mafdr(post_inj,'BH',1);
pfdr = [ones(1,238) pfdr];
idx = double(find(pfdr<0.05));
hold on
    plot(mean(fpn_CE)','black','LineWidth',2)
    plot(mean(fpn_CE_pcb)','LineWidth',2)
text(idx',repelem(.9*max(mean(fpn_CE_pcb)),length(idx)),'*')
title(networks{6})
% legend('DMT','PCB')
xticks(tics)
xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
xlabel('Minutes')


subplot(5,2,9)
[h,p]=ttest(dmn_CE,dmn_CE_pcb);
% idx = double(find(p<0.05));
post_inj = p(239:end);
pfdr = mafdr(post_inj,'BH',1);
pfdr = [ones(1,238) pfdr];
idx = double(find(pfdr<0.05));
hold on
    plot(mean(dmn_CE)','black','LineWidth',2)
    plot(mean(dmn_CE_pcb)','LineWidth',2)
text(idx',repelem(.9*max(mean(dmn_CE_pcb)),length(idx)),'*')
title(networks{7})
% legend('DMT','PCB')
xticks(tics)
xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
xlabel('Minutes')


subplot(5,2,10)
[h,p]=ttest(sub_CE,sub_CE_pcb);
% idx = double(find(p<0.05));
post_inj = p(239:end);
pfdr = mafdr(post_inj,'BH',1);
pfdr = [ones(1,238) pfdr];
idx = double(find(pfdr<0.05));
hold on
    plot(mean(sub_CE)','black','LineWidth',2)
    plot(mean(sub_CE_pcb)','LineWidth',2)
text(idx',repelem(.9*max(mean(sub_CE_pcb)),length(idx)),'*')
title('SUB')
% legend('DMT','PCB')
xticks(tics)
xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
xlabel('Minutes')


%% SI FIGURE 2 summarize network comparison
[h,p,~,t]=ttest(global_CE_dmt,global_CE_pcb);
post_inj_glob = t.tstat(239:end);
post_inj1_glob = post_inj_glob(1:300);
post_inj2_glob = post_inj_glob(301:600);

[h,p,~,t]=ttest(vis_CE,vis_CE_pcb);
post_inj_vis = t.tstat(239:end);
post_inj1_vis = post_inj_vis(1:300);
post_inj2_vis = post_inj_vis(301:600);

[h,p,~,t]=ttest(som_CE,som_CE_pcb);
post_inj_som = t.tstat(239:end);
post_inj1_som = post_inj_som(1:300);
post_inj2_som = post_inj_som(301:600);

[h,p,~,t]=ttest(dat_CE,dat_CE_pcb);
post_inj_dat = t.tstat(239:end);
post_inj1_dat = post_inj_dat(1:300);
post_inj2_dat = post_inj_dat(301:600);

[h,p,~,t]=ttest(vat_CE,vat_CE_pcb);
post_inj_vat = t.tstat(239:end);
post_inj1_vat = post_inj_vat(1:300);
post_inj2_vat = post_inj_vat(301:600);

[h,p,~,t]=ttest(lim_CE,lim_CE_pcb);
post_inj_lim = t.tstat(239:end);
post_inj1_lim = post_inj_lim(1:300);
post_inj2_lim = post_inj_lim(301:600);

[h,p,~,t]=ttest(fpn_CE,fpn_CE_pcb);
post_inj_fpn = t.tstat(239:end);
post_inj1_fpn = post_inj_fpn(1:300);
post_inj2_fpn = post_inj_fpn(301:600);

[h,p,~,t]=ttest(dmn_CE,dmn_CE_pcb);
post_inj_dmn = t.tstat(239:end);
post_inj1_dmn = post_inj_dmn(1:300);
post_inj2_dmn = post_inj_dmn(301:600);

[h,p,~,t]=ttest(dmn_CE,dmn_CE_pcb);
post_inj_dmn = t.tstat(239:end);
post_inj1_dmn = post_inj_dmn(1:300);
post_inj2_dmn = post_inj_dmn(301:600);

[h,p,~,t]=ttest(sub_CE,sub_CE_pcb);
post_inj_sub = t.tstat(239:end);
post_inj1_sub = post_inj_sub(1:300);
post_inj2_sub = post_inj_sub(301:600);

% data = [post_inj1_glob' post_inj2_glob' post_inj1_vis' post_inj2_vis' post_inj1_som' post_inj2_som' post_inj1_dat' post_inj2_dat' post_inj1_vat' post_inj2_vat' post_inj1_lim' post_inj2_lim' post_inj1_fpn' post_inj2_fpn' post_inj1_dmn' post_inj2_dmn' post_inj1_sub' post_inj2_sub'];
data = [post_inj1_vis' post_inj2_vis' post_inj1_fpn' post_inj2_fpn' post_inj1_dmn' post_inj2_dmn'];
figure;
violin(data)
xticks(1:6);
xticklabels([{'VIS first half'},{'VIS second half'},{'FPN first half'},{'FPN second half'},{'DMN first half'},{'DMN second half'}]);
xtickangle(45)
ylabel('t-stat (DMT > PCB)')



%% FIGURE 2a
figure;
[h,p]=ttest(global_CE_dmt,global_CE_pcb);
post_inj = p(239:end);
pfdr = mafdr(post_inj,'BH',1);
fraction_sig_after_injection = sum(pfdr<0.05)/length(post_inj)
pfdr = [ones(1,238) pfdr];
idx = double(find(pfdr<0.05));
hold on
    plot(mean(global_CE_dmt)','black','LineWidth',1.5)
    plot(mean(global_CE_pcb)','red','LineWidth',1.5)
text(idx',repelem(.9*max(mean(global_CE_pcb)),length(idx)),'*')
tics = linspace(0,28,15);
tics = tics*30;
tics(end)=838;
xticks(tics)
xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
xlabel('Minutes')
ylabel('Control Energy');
title([{'Group-level continuous control energy'};{'DMT vs PCB'}])
legend('DMT Control Energy','PCB Control Energy')

%% window global CE for intensity comparison

load([basedir,'data/intensity_ratings.mat'])

m_dmt_int = mean(dmt_intensity);
m_pcb_int = mean(pcb_intensity);


TR=2; window=60/TR;

win_dmt_ce=[];
win_pcb_ce=[];
win_dmt_ce_global=[];
win_pcb_ce_global=[];

%get avg energy for each min
k=1; 
for i=1:length(m_dmt_int)
    if i==length(m_dmt_int)
        win_dmt_ce(:,i,:) = squeeze(mean(regional_CE_dmt(:,k:end,:),2));
        win_pcb_ce(:,i,:) = squeeze(mean(regional_CE_pcb(:,k:end,:),2));
        win_dmt_ce_global(:,i) = squeeze(mean(global_CE_dmt(:,k:end),2));
        win_pcb_ce_global(:,i) = squeeze(nanmean(global_CE_pcb(:,k:end),2));

    else
        win_dmt_ce(:,i,:) = squeeze(mean(regional_CE_dmt(:,k:k+window,:),2));
        win_pcb_ce(:,i,:) = squeeze(mean(regional_CE_pcb(:,k:k+window,:),2));
        win_dmt_ce_global(:,i,:) = squeeze(mean(global_CE_dmt(:,k:k+window),2));
        win_pcb_ce_global(:,i,:) = squeeze(mean(global_CE_pcb(:,k:k+window),2));

    end
    k=k+window;
end


%% CE vs intensity for each condition

nperms=10000;

clear rr pp

figure;
subplot(2,1,1)
data = mean(win_dmt_ce_global);
[r,p] = corr(data',m_dmt_int','type','Spearman');
for i=1:nperms
    idx1 = randperm(length(data));
%     idx2 = randperm(length(data));
    rr(i) = corr(data(idx1)',m_dmt_int','type','Spearman');
end
pp = mean(r>rr);
title([{'Energy/Intensity Time-series'};{'DMT continuous data (group means)'}])
hold on
    plot(data,'bo--','LineWidth',2)
    text(14,160,['R = ',num2str(r),'; p = ',num2str(pp)],'FontSize',12)
    ylabel('E')
%     ylim([250 950])
    ylim([50 180])
    yyaxis right
    plot(m_dmt_int,'blacko--','LineWidth',2)
ylabel('Intensity')
ylim([0 10])
xlabel('Minutes')
xticks(linspace(0,28,15))
xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
legend('Energy','Intensity')

subplot(2,1,2)
data=mean(win_pcb_ce_global);
[r,p]=corr(data',m_pcb_int','type','Spearman')
for i=1:nperms
    idx1 = randperm(length(data));
%     idx2 = randperm(length(data));
    rr(i) = corr(data(idx1)',m_pcb_int','type','Spearman');
end
pp = mean(r>rr);
title([{'Energy/Intensity Time-series'};{'PCB continuous data (group means)'}])
hold on
    plot(data,'ro--','LineWidth',2)
    text(14,160,['R = ',num2str(r),'; p = ',num2str(pp)],'FontSize',12)
    ylabel('E')
%     ylim([250 950])
    ylim([50 180])
    yyaxis right
    plot(m_pcb_int,'blacko--','LineWidth',2)
ylabel('Intensity')
ylim([0 10])
xlabel('Minutes')
xticks(linspace(0,28,15))
xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
legend('Energy','Intensity')

%% FIGURE 2B and 2C

load([basedir,'data/RegressorLZInterpscrubbedConvolvedAvg.mat'])

%baseline correct EEG first:
bldmt = nanmean(RegDMT2(:,1:240),2);
RegDMT2 = RegDMT2 - bldmt;
RegPCB2 = RegPCB2 - nanmean(RegPCB2(:,1:240),2);

RegDMT2 = RegDMT2(:,2:839);
m_dmt_LZ = nanmean(RegDMT2);
RegPCB2 = RegPCB2(:,2:839);
m_pcb_LZ = nanmean(RegPCB2);

diff_ce = global_CE_dmt - global_CE_pcb;
diff_lz = RegDMT2 - RegPCB2;

m_diff_lz = mean(diff_lz);

figure;
subplot(2,1,1)
hold on
    plot(nanmean(diff_ce)','blue','LineWidth',1.5)
    
    ylabel('∆Control Energy')
    yyaxis right
    plot(mean(diff_lz)','green','LineWidth',1.5)
    ylabel('∆EEG Signal Diversity (LZ)')
    yyaxis left
[r,p] = corr(nanmean(diff_ce)',m_diff_lz','type','Spearman')
for i=1:nperms
    idx1 = randperm(length(m_dmt_LZ));
    rr(i) = corr(nanmean(diff_ce)',m_diff_lz(idx1)','type','Spearman');
end
pp = mean(r>=rr)
text(400,-100,['R = ',num2str(r),'; p = ',num2str(pp)],'FontSize',12)
title([{'Group-level continuous control energy vs signal diversity'};{'(DMT - PCB)'}])
xlabel('Minutes')
tics = linspace(0,28,15);
tics = tics*30;
tics(end)=838;
xticks(tics)
xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
legend('∆Control Energy (fMRI)','∆Limpel-Ziv Complexity (EEG)')

clear rr pp

subplot(2,1,2)
data = mean(win_dmt_ce_global-win_pcb_ce_global);
diff_int = m_dmt_int - m_pcb_int;
[r,p] = corr(data',diff_int','type','Spearman');
for i=1:nperms
    idx1 = randperm(length(data));
    rr(i) = corr(data(idx1)',diff_int','type','Spearman');
end
pp = mean(r>=rr)
title([{'Group-level windowed control energy vs drug intensity'};{'(DMT - PCB)'}])
hold on
    plot(data,'bo--','LineWidth',1.5)
    text(18,-30,['R = ',num2str(r),'; p = ',num2str(pp)],'FontSize',12)
    ylabel('∆Control Energy')
    yyaxis right
    plot(diff_int,'yellowo--','LineWidth',1.5)
ylabel('∆Intensity Ratings')
ylim([0 10])
xlabel('Minutes')
xticks(linspace(0,28,15))
xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
legend('∆Control Energy (fMRI)','∆Intensity (subjective ratings)')

