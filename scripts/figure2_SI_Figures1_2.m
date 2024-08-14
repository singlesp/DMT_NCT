%% plot global and network level CE - must run <gen_time_resolved_ce.m> first

clear all; close all;

nsub=14;

nperms=10000;

basedir = '~/Documents/GIT/DMT_NCT/';

note='_gsr_volnorm'; %track different proc streams

load([basedir,'results/regional_continuous_CE_DMT',note,'.mat'])
load([basedir,'results/regional_continuous_CE_PCB',note,'.mat'])

load([basedir,'data/sch116_to_yeo.csv'])
network7labels=sch116_to_yeo;
networks = [{'VIS'},{'SOM'},{'DAT'},{'VAT'},{'LIM'},{'FPN'},{'DMN'},{'SUB'}];

%% get network level CE's for SI

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

%% SI FIGURE Network level

% tics = linspace(0,28,15);
% tics = tics*30;
% tics(end)=838;


figure;
subplot(5,2,1.5)
run_cluster_stats(global_CE_dmt,global_CE_pcb,-1,{'CE'},[1 0 0],[0 0 0],[],[],0,{'Global'},[12 14 18]);

subplot(5,2,3)
run_cluster_stats(vis_CE,vis_CE_pcb,-1,{'CE'},[1 0 0],[0 0 0],[],[],0,networks{1},[12 14 18]);

subplot(5,2,4)
run_cluster_stats(som_CE,som_CE_pcb,-1,{'CE'},[1 0 0],[0 0 0],[],[],0,networks{2},[12 14 18]);

subplot(5,2,5)
run_cluster_stats(dat_CE,dat_CE_pcb,-1,{'CE'},[1 0 0],[0 0 0],[],[],0,networks{3},[12 14 18]);

subplot(5,2,6)
run_cluster_stats(vat_CE,vat_CE_pcb,-1,{'CE'},[1 0 0],[0 0 0],[],[],0,networks{4},[12 14 18]);

subplot(5,2,7)
run_cluster_stats(lim_CE,lim_CE_pcb,-1,{'CE'},[1 0 0],[0 0 0],[],[],0,networks{5},[12 14 18]);

subplot(5,2,8)
run_cluster_stats(fpn_CE,fpn_CE_pcb,-1,{'CE'},[1 0 0],[0 0 0],[],[],0,networks{6},[12 14 18]);

subplot(5,2,9)
run_cluster_stats(dmn_CE,dmn_CE_pcb,-1,{'CE'},[1 0 0],[0 0 0],[],[],0,networks{7},[12 14 18]);

subplot(5,2,10)
run_cluster_stats(sub_CE,sub_CE_pcb,-1,{'CE'},[1 0 0],[0 0 0],[],[],0,{'SUB'},[12 14 18]);

%% SI FIGURE summarize network comparison
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
run_cluster_stats(global_CE_dmt,global_CE_pcb,-1,{'CE'},[1 0 0],[0 0 0]);


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


%% CE vs intensity for each condition (SI)

clear rr pp*

figure;
subplot(2,1,1)
data = mean(win_dmt_ce_global);
[r,p] = corr(data',m_dmt_int','type','Spearman');
for i=1:nperms
    idx1 = randperm(length(data));
    rr(i) = corr(data(idx1)',m_dmt_int','type','Spearman');
end
pp1_int = mean(r>rr);
title([{'Energy/Intensity Time-series'};{'DMT continuous data (group means)'}])
hold on
    plot(data,'bo--','LineWidth',2)
    text(14,160,['R = ',num2str(r),'; p = ',num2str(pp1_int)],'FontSize',12)
    ylabel('Global Control Energy (a.u.)')
%     ylim([250 950])
    ylim([50 180])
    yyaxis right
    plot(m_dmt_int,'blacko--','LineWidth',2)
ylabel('Intensity Ratings')
ylim([0 10])
xlabel('Minutes')
xlim([0 28])
xticks(linspace(0,28,8))
xticklabels([{'-8'},{'-4'},{'0'},{'4'},{'8'},{'12'},{'16'},{'20'}]);
legend('Energy','Intensity')

subplot(2,1,2)
data=mean(win_pcb_ce_global);
[r,p]=corr(data',m_pcb_int','type','Spearman');
for i=1:nperms
    idx1 = randperm(length(data));
    rr(i) = corr(data(idx1)',m_pcb_int','type','Spearman');
end
pp2_int = mean(r>rr);
title([{'Energy/Intensity Time-series'};{'PCB continuous data (group means)'}])
hold on
    plot(data,'ro--','LineWidth',2)
    text(14,160,['R = ',num2str(r),'; p = ',num2str(pp2_int)],'FontSize',12)
    ylabel('Global Control Energy (a.u.)')
%     ylim([250 950])
    ylim([50 180])
    yyaxis right
    plot(m_pcb_int,'blacko--','LineWidth',2)
ylabel('Intensity Ratings')
ylim([0 10])
xlabel('Minutes')
xlim([0 28])
xticks(linspace(0,28,8))
xticklabels([{'-8'},{'-4'},{'0'},{'4'},{'8'},{'12'},{'16'},{'20'}]);
legend('Energy','Intensity')

%% CE vs LZ for each condition

load([basedir,'data/RegressorLZInterpscrubbedConvolvedAvg.mat'])

%baseline correct EEG first:
bldmt = nanmean(RegDMT2(:,1:240),2);
RegDMT2 = RegDMT2 - bldmt;
RegPCB2 = RegPCB2 - nanmean(RegPCB2(:,1:240),2);

RegDMT2 = RegDMT2(:,2:839);
m_dmt_LZ = nanmean(RegDMT2);
RegPCB2 = RegPCB2(:,2:839);
m_pcb_LZ = nanmean(RegPCB2);

m_dmt_ce = nanmean(global_CE_dmt);
m_pcb_ce = nanmean(global_CE_pcb);


figure;
subplot(2,1,1)
hold on
    plot(m_dmt_ce','blue','LineWidth',1.5)
    
    ylabel('Global Control Energy (a.u.)')
    yyaxis right
    plot(m_dmt_LZ','green','LineWidth',1.5)
    ylabel('EEG Signal Diversity (LZ)')
    yyaxis left
[r,p] = corr(m_dmt_ce',m_dmt_LZ','type','Spearman');
for i=1:nperms
    idx1 = randperm(length(m_dmt_LZ));
    rr(i) = corr(m_dmt_ce',m_dmt_LZ(idx1)','type','Spearman');
end
pp1_eg = mean(r>=rr)
text(350,100,['R = ',num2str(r),'; p = ',num2str(pp1_eg)],'FontSize',12)
title([{'DMT continuous control energy vs signal diversity'}])
xlabel('Minutes')
xlim([0 838])
tics = linspace(0,28,8);
tics = tics*30;
tics(end)=838;
xticks(tics)
xticklabels([{'-8'},{'-4'},{'0'},{'4'},{'8'},{'12'},{'16'},{'20'}]);
legend('Control Energy (fMRI)','Limpel-Ziv Complexity (EEG)')


subplot(2,1,2)
hold on
    plot(m_pcb_ce','blue','LineWidth',1.5)
    
    ylabel('Global Control Energy (a.u.)')
    yyaxis right
    plot(m_pcb_LZ','green','LineWidth',1.5)
    ylabel('EEG Signal Diversity (LZ)')
    yyaxis left
[r,p] = corr(m_pcb_ce',m_pcb_LZ','type','Spearman');
for i=1:nperms
    idx1 = randperm(length(m_pcb_LZ));
    rr(i) = corr(m_pcb_ce',m_pcb_LZ(idx1)','type','Spearman');
end
pp2_eg = mean(r>=rr)
text(300,200,['R = ',num2str(r),'; p = ',num2str(pp2_eg)],'FontSize',12)
title([{'PCB continuous control energy vs signal diversity'}])
xlabel('Minutes')
xlim([0 838])
tics = linspace(0,28,8);
tics = tics*30;
tics(end)=838;
xticks(tics)
xticklabels([{'-8'},{'-4'},{'0'},{'4'},{'8'},{'12'},{'16'},{'20'}]);
legend('Control Energy (fMRI)','Limpel-Ziv Complexity (EEG)')

%% correct global single condition SI correlations

pfdr = mafdr([pp1_int pp2_int pp1_eg pp2_eg],'BH',1)


%% FIGURE 2B and 2C

clear rr pp*

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
pp1 = mean(r>=rr)
text(400,-100,['R = ',num2str(r),'; p = ',num2str(pp1)],'FontSize',12)
title([{'Group-level continuous control energy vs signal diversity'};{'(DMT - PCB)'}])
xlabel('Minutes')
tics = linspace(0,28,15);
tics = tics*30;
tics(end)=838;
xlim([0 838])
xticks(tics)
xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
legend('∆Control Energy (fMRI)','∆Limpel-Ziv Complexity (EEG)')


subplot(2,1,2)
data = mean(win_dmt_ce_global-win_pcb_ce_global);
diff_int = m_dmt_int - m_pcb_int;
[r,p] = corr(data',diff_int','type','Spearman');
for i=1:nperms
    idx1 = randperm(length(data));
    rr(i) = corr(data(idx1)',diff_int','type','Spearman');
end
pp2 = mean(r>=rr)
title([{'Group-level windowed control energy vs drug intensity'};{'(DMT - PCB)'}])
hold on
    plot(data,'bo--','LineWidth',1.5)
    text(18,-30,['R = ',num2str(r),'; p = ',num2str(pp2)],'FontSize',12)
    ylabel('∆Control Energy')
    yyaxis right
    plot(diff_int,'yellowo--','LineWidth',1.5)
ylabel('∆Intensity Ratings')
ylim([0 10])
xlabel('Minutes')
xlim([0 28])
xticks(linspace(0,28,15))
xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
legend('∆Control Energy (fMRI)','∆Intensity (subjective ratings)')


%% correlation for post-injection period only (in-text)

post_inj_diff_lz = m_diff_lz(239:end);
[r,p] = corr(nanmean(diff_ce(:,239:end))',post_inj_diff_lz','type','Spearman')
for i=1:nperms
    idx1 = randperm(length(post_inj_diff_lz));
    rr(i) = corr(nanmean(diff_ce(:,239:end))',post_inj_diff_lz(idx1)','type','Spearman');
end
pp1_post = mean(r>=rr);


data = mean(win_dmt_ce_global-win_pcb_ce_global);
data_post = data(9:end);
diff_int = m_dmt_int - m_pcb_int;
diff_int_post = diff_int(9:end);
[r,p] = corr(data_post',diff_int_post','type','Spearman')
for i=1:nperms
    idx1 = randperm(length(data_post));
    rr(i) = corr(data_post(idx1)',diff_int_post','type','Spearman');
end
pp2_post = mean(r>=rr);

%% in-text correlation for full scan: partial corr w/ FD covar

load([basedir,'data/FDlong.mat'])

m_dmt_fd = mean(FDDMT(2:839,:)');
m_pcb_fd = mean(FDPCB(2:839,:)');
diff_fd = m_dmt_fd - m_pcb_fd;

% window the fd for correlation with intensity
TR=2; window=60/TR;

win_dmt_fd=[];
win_pcb_fd=[];

%get avg energy for each min
k=1; 
for i=1:length(m_dmt_int)
    if i==length(m_dmt_int)
        win_dmt_fd(:,i) = mean(FDDMT(k:end,:)',2);
        win_pcb_fd(:,i) = mean(FDPCB(k:end,:)',2);
    else
        win_dmt_fd(:,i) = mean(FDDMT(k:k+window,:)',2);
        win_pcb_fd(:,i) = mean(FDPCB(k:k+window,:)',2);

    end
    k=k+window;
end

% partial spearman corr bw energy diff and LZ
resid = regress_confound([tiedrank(nanmean(diff_ce)') tiedrank(m_diff_lz')],tiedrank(diff_fd'),'fittype','lsq','addconstant',true);
[r,p] = corr(resid(:,1),resid(:,2))
for i=1:nperms
    idx1 = randperm(length(m_dmt_LZ));
    rr(i) = corr(resid(:,1),resid(idx1,2));
end
pp1_fd = mean(r>=rr)

% partial spearman for bw energy diff and intensity 
m_win_ce_diff = mean(win_dmt_ce_global-win_pcb_ce_global);
m_win_fd_diff = mean(win_dmt_fd-win_pcb_fd);
diff_int = m_dmt_int - m_pcb_int;

resid = regress_confound([tiedrank(m_win_ce_diff') tiedrank(diff_int')],tiedrank(m_win_fd_diff'),'fittype','lsq','addconstant',true);
[r,p] = corr(resid(:,1),resid(:,2))
for i=1:nperms
    idx1 = randperm(length(diff_int));
    rr(i) = corr(resid(:,1),resid(idx1,2));
end
pp2_fd = mean(r>=rr)

%% correction of main text global correlations

pfdr = mafdr([pp1 pp2 pp1_post pp2_post pp1_fd pp2_fd],'BH',1)

%% EEG direct comparison (SI)

load([basedir,'data/RegressorLZInterpscrubbedConvolvedAvg.mat'])

%baseline correct EEG first:
bldmt = nanmean(RegDMT2(:,1:240),2);
RegDMT2 = RegDMT2 - bldmt;
RegPCB2 = RegPCB2 - nanmean(RegPCB2(:,1:240),2);

RegDMT2 = RegDMT2(:,2:839);
m_dmt_LZ = nanmean(RegDMT2);
RegPCB2 = RegPCB2(:,2:839);
m_pcb_LZ = nanmean(RegPCB2);

% diff_ce = global_CE_dmt - global_CE_pcb;
% diff_lz = RegDMT2 - RegPCB2;

figure;
[h,p]=ttest(RegDMT2,RegPCB2);
post_inj = p(239:end);
pfdr = mafdr(post_inj,'BH',1);
fraction_sig_after_injection = sum(pfdr<0.05)/length(post_inj)
pfdr = [ones(1,238) pfdr];
idx = double(find(pfdr<0.05));
hold on
    plot(m_dmt_LZ,'black','LineWidth',1.5)
    plot(m_pcb_LZ,'red','LineWidth',1.5)
text(idx',repelem(.9*max(m_dmt_LZ),length(idx)),'*')
tics = linspace(0,28,15);
tics = tics*30;
tics(end)=838;
xticks(tics)
xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
xlabel('Minutes')
ylabel('EEG Signal Diversity (LZ)');
title([{'Entropy time-series'};{'DMT vs PCB'}])
legend('DMT signal diversity','PCB signal diversity')

%% subject level correlations (SI)
%DMT CE vs LZ

for i=1:nsub
    [r,p] = corr(global_CE_dmt(i,:)',RegDMT2(i,:)','type','Spearman');
    for perm=1:nperms
        idx1 = randperm(length(diff_lz));
        rr(perm) = corr(global_CE_dmt(i,:)',RegDMT2(i,idx1)','type','Spearman');
    end
    if r<0
        pp(i) = mean(r>=rr);
    else
        pp(i) = mean(r<=rr);
    end
end
pfdr = mafdr(pp,'BH',1);

figure;
for i=1:nsub
    subplot(7,2,i)
    hold on
    plot(global_CE_dmt(i,:)','blue','LineWidth',1.5)
    
    ylabel('CE')
    yyaxis right
    plot(RegDMT2(i,:)','green','LineWidth',1.5)
    ylabel('LZ')
    yyaxis left
    [r,p] = corr(global_CE_dmt(i,:)',RegDMT2(i,:)','type','Spearman');
    text(400,200,['R = ',num2str(round(r,2)),'; p = ',num2str(round(pfdr(i),4))],'FontSize',12)
    xlabel('Minutes')
    tics = linspace(0,28,15);
    tics = tics*30;
    tics(end)=838;
    xlim([0 838])
    xticks(tics)
    xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
end
set(gcf,'Position',[0 0 1000 975])


%% PCB CE vs LZ

for i=1:nsub
   [r,p] = corr(global_CE_pcb(i,:)',RegPCB2(i,:)','type','Spearman','rows','complete');
    for perm=1:nperms
        idx1 = randperm(length(diff_lz));
        rr(perm) = corr(global_CE_pcb(i,:)',RegPCB2(i,idx1)','type','Spearman','rows','complete');
    end
    if r<0
        pp(i) = mean(r>=rr);
    else
        pp(i) = mean(r<=rr);
    end 
end
pfdr = mafdr(pp,'BH',1);

figure;
for i=1:nsub
    subplot(7,2,i)
    hold on
    plot(global_CE_pcb(i,:)','blue','LineWidth',1.5)
    
    ylabel('CE')
    yyaxis right
    plot(RegPCB2(i,:)','green','LineWidth',1.5)
    ylabel('LZ')
    yyaxis left
    [r,p] = corr(global_CE_pcb(i,:)',RegPCB2(i,:)','type','Spearman','rows','complete');
    text(400,200,['R = ',num2str(round(r,2)),'; p = ',num2str(round(pfdr(i),4))],'FontSize',12)
    xlabel('Minutes')
    tics = linspace(0,28,15);
    tics = tics*30;
    tics(end)=838;
    xlim([0 838])
    xticks(tics)
    xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
end
set(gcf,'Position',[0 0 1000 975])

%% DMT CE vs intensity

for i=1:nsub
    [r,p] = corr(win_dmt_ce_global(i,:)',dmt_intensity(i,:)','type','Spearman');
   for perm=1:nperms
        idx1 = randperm(length(dmt_intensity));
        rr(perm) = corr(win_dmt_ce_global(i,idx1)',dmt_intensity(i,:)','type','Spearman');
    end
    if r<0
        pp2(i) = mean(r>=rr);
    else
        pp2(i) = mean(r<=rr);
    end 
end

pfdr = mafdr(pp2,'BH',1);

figure;
for i=1:nsub
    subplot(7,2,i)
    [r,p] = corr(win_dmt_ce_global(i,:)',dmt_intensity(i,:)','type','Spearman');
    hold on
    plot(win_dmt_ce_global(i,:)','bo--','LineWidth',1.5)
    text(16,100,['R = ',num2str(round(r,2)),'; p = ',num2str(round(pfdr(i),4))],'FontSize',12)
    ylabel('CE')
    yyaxis right
    plot(dmt_intensity(i,:)','go--','LineWidth',1.5)
    ylabel('Intensity Ratings')
    ylim([0 10])
    xlabel('Minutes')
    xlim([0 28])
    xticks(linspace(0,28,15))
    xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
end
set(gcf,'Position',[0 0 1000 975])

%%

for i=1:nsub
   [r,p] = corr(win_pcb_ce_global(i,:)',pcb_intensity(i,:)','type','Spearman');
    for perm=1:nperms
        idx1 = randperm(length(pcb_intensity));
        rr(perm) = corr(win_pcb_ce_global(i,idx1)',pcb_intensity(i,:)','type','Spearman');
    end
    if r<0
        pp2(i) = mean(r>=rr);
    else
        pp2(i) = mean(r<=rr);
    end  
end
pfdr = mafdr(pp2,'BH',1);

figure;
for i=1:nsub
    subplot(7,2,i)
    [r,p] = corr(win_pcb_ce_global(i,:)',pcb_intensity(i,:)','type','Spearman');
    hold on
    plot(win_pcb_ce_global(i,:)','bo--','LineWidth',1.5)
    text(16,100,['R = ',num2str(round(r,2)),'; p = ',num2str(round(pfdr(i),4))],'FontSize',12)
    ylabel('CE')
    yyaxis right
    plot(pcb_intensity(i,:)','go--','LineWidth',1.5)
    ylabel('Intensity Ratings')
    ylim([0 10])
    xlabel('Minutes')
    xlim([0 28])
    xticks(linspace(0,28,15))
    xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
end
set(gcf,'Position',[0 0 1000 975])

%%


figure;
for i=1:nsub
    subplot(7,2,i)
    cei = global_CE_pcb(i,:)-global_CE_dmt(i,:);
    sdi = RegPCB2(i,:)-RegDMT2(i,:);
    hold on
    plot(cei','blue','LineWidth',1.5)
    
    ylabel('CE')
    yyaxis right
    plot(sdi','green','LineWidth',1.5)
    ylabel('LZ')
    yyaxis left
    [r,p] = corr(cei',sdi','type','Spearman');
    for perm=1:nperms
        idx1 = randperm(length(diff_lz));
        rr(perm) = corr(cei',sdi(1,idx1)','type','Spearman');
    end
    pp(i) = mean(r>=rr)
    text(400,200,['R = ',num2str(r),'; p = ',num2str(pp(i))],'FontSize',12)
%     title([{'Group-level continuous control energy vs signal diversity'};{'(DMT - PCB)'}])
    xlabel('Minutes')
    tics = linspace(0,28,15);
    tics = tics*30;
    tics(end)=838;
    xlim([0 838])
    xticks(tics)
    xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
%     legend('∆Control Energy (fMRI)','∆Limpel-Ziv Complexity (EEG)')
end

pfdr = mafdr(pp,'BH',1)
