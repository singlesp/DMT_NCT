%% regional CE for mins8 (Figure 3, left)
clear all; close all;

basedir = '~/Documents/GIT/DMT_NCT/';

note='_gsr_volnorm'; %track different proc streams

nsub=14;


%% corr with sub included (not reported in paper)

load(fullfile([basedir,'results/regional_continuous_prepost_CE_PCB',note,'.mat']),'global_CE_pcb','regional_CE_pcb','global_CE_pcb_pre','regional_CE_pcb_pre');

load([basedir,'data/5HTvecs_sch116.mat'], 'mean5HT2A_sch116');
receptor_vec = zscore(mean5HT2A_sch116);

pre_energy = squeeze(mean(regional_CE_pcb_pre,2));
post_energy = squeeze(mean(regional_CE_pcb,2));

diff_energy = (post_energy-pre_energy)./(pre_energy);

[r2a,p2a] = corr(mean(diff_energy)',receptor_vec,'type','Spearman')


figure;
scatter(receptor_vec,mean(diff_energy)'); lsline
xlabel('2a density')
ylabel('pre/post PCB Nodal Control Energy (∆)')
text(-2,-0.6,['rho = ',num2str(r2a),'; p = ',num2str(p2a)])

%% FIGURE 3a, left.  gummibrain from Keith Jamison-> https://github.com/kjamison/atlasblobs

load(fullfile([basedir,'fxns/gummi_kj/atlasblobs_saved.mat']))

%choose atlas
whichatlas={'sch116'};

%set colormap
cmap = brewermap([],'*RdBu');

data = mean(diff_energy)';
clim=[-0.6 0.6];
 
f=figure;

img=display_atlas_blobs(data,atlasblobs,...
    'atlasname',whichatlas,...
    'render',true,...
    'backgroundimage',true,...
    'crop',true,...
    'colormap',cmap,...
    'clim', clim,...
    'backgroundimage','none',...
    'backgroundcolor',[1 1 1]);%,...
%         'roimask',roimask); %last argmument optional, depends if you need to remove roi's

imshow(img);
c=colorbar('SouthOutside', 'fontsize', 16);
c.Label.String='∆CE: (post - pre)/pre';
set(gca,'colormap',cmap);
caxis(clim);

%% FIGURE 3c, left (same corr as above but no subcortex so can spin)

load([basedir,'fxns/SpinTests/rotated_maps/rotated_Schaefer_100.mat'])

diff_energy= diff_energy(:,1:100);
rec_nosub = receptor_vec(1:100);

[r2a,p2a] = corr(mean(diff_energy)',rec_nosub,'type','Spearman')

p_spin_vasa_left = perm_sphere_p_al857(mean(diff_energy)',rec_nosub,perm_id,'Spearman')

figure;
scatter(rec_nosub,mean(diff_energy)'); lsline
xlabel('2a density')
ylabel('pre/post PCB Nodal Control Energy (∆)')
text(-.7,0.15,['rho = ',num2str(r2a),'; p = ',num2str(p_spin_vasa_left)])


%% dominance analysis - see figure4.m for plotting these data
load([basedir,'data/5HTvecs_sch116.mat'])

pre_energy = squeeze(mean(regional_CE_pcb_pre,2));
post_energy = squeeze(mean(regional_CE_pcb,2));

diff_energy = (post_energy-pre_energy)./(pre_energy);

v2a=mean5HT2A_sch116;

v1a=mean5HT1A_sch116;

v1b=mean5HT1B_sch116;

v4=mean5HT4_sch116;

vT=mean5HTT_sch116;

allv = [v2a v1a v1b v4 vT];

target = mean(diff_energy)';

varnames=[{'HT2a'},{'HT1a'},{'HT1b'},{'HT4'},{'HTT'},{'target'}];

mydata = table(v2a,v1a,v1b,v4,vT,target,'VariableNames',varnames);

input_filename = 'results/dominance/diffE_pcb_continuous';
writetable(mydata, [basedir,input_filename,note, '.csv']);

%call from terminal:
% python3 ~/Documents/GIT/DMT_NCT/fxns/dominance_analysis_al857.py ~/Documents/GIT/DMT_NCT/results/dominance/diffE_pcb_continuous_gsr_volnorm.csv  target 5 0 ~/Documents/GIT/DMT_NCT/results/dominance/diffE_pcb_continuous_gsr_volnorm

%% regional CE vs LZ; corr w sub

load([basedir,'data/RegressorLZInterpscrubbedConvolvedAvg.mat'])

m_pcb_LZ = nanmean(RegPCB2);
m_pcb_LZ = m_pcb_LZ(2:839); %first frame is excluded from CE calcs


[r_avg,p_avg]=corr(m_pcb_LZ',squeeze(nanmean(regional_CE_pcb,1)),'type','Spearman');


[r2a,p2a] = corr(r_avg',receptor_vec,'type','Spearman')
figure;
scatter(receptor_vec,r_avg'); lsline
xlabel('2a density')
ylabel('Rho: Nodal correlation b/w CE and LZ')
text(-2,0.2,['rho = ',num2str(r2a),'; p = ',num2str(p2a)])

%% FIGURE 3a, middle.  gummibrain from Keith Jamison-> https://github.com/kjamison/atlasblobs
load(fullfile([basedir,'fxns/gummi_kj/atlasblobs_saved.mat']))

%choose atlas
whichatlas={'sch116'};

%set colormap
cmap = brewermap([],'*RdBu');

data = r_avg';
clim=[-0.6 0.6];
 
f=figure;

img=display_atlas_blobs(data,atlasblobs,...
    'atlasname',whichatlas,...
    'render',true,...
    'backgroundimage',true,...
    'crop',true,...
    'colormap',cmap,...
    'clim', clim,...
    'backgroundimage','none',...
    'backgroundcolor',[1 1 1]);%,...
%         'roimask',roimask); %last argmument optional, depends if you need to remove roi's

imshow(img);
c=colorbar('SouthOutside', 'fontsize', 16);
c.Label.String='∆CE: (post - pre)/pre';
set(gca,'colormap',cmap);
caxis(clim);

%% FIGURE 3c, middle above correlation without subcortex and using spin test
load([basedir,'fxns/SpinTests/rotated_maps/rotated_Schaefer_100.mat'])


r_avg=r_avg(1:100);
rec_nosub = receptor_vec(1:100);

[r2a,p2a] = corr(r_avg',rec_nosub,'type','Spearman')

p_spin_vasa_mid = perm_sphere_p_al857(r_avg',rec_nosub,perm_id,'Spearman')

figure;
scatter(rec_nosub,r_avg'); lsline
xlabel('2a density')
ylabel('Rho: Nodal correlation b/w CE and LZ')
text(-.7,0.1,['rho = ',num2str(r2a),'; p = ',num2str(p_spin_vasa_mid)])

%% dominance analysis - see figure4.m for plotting these data
load([basedir,'data/5HTvecs_sch116.mat'])

v2a=mean5HT2A_sch116;

v1a=mean5HT1A_sch116;

v1b=mean5HT1B_sch116;

v4=mean5HT4_sch116;

vT=mean5HTT_sch116;

allv = [v2a v1a v1b v4 vT];

[r_avg,p_avg]=corr(m_pcb_LZ',squeeze(nanmean(regional_CE_pcb,1)));
target = r_avg';

varnames=[{'HT2a'},{'HT1a'},{'HT1b'},{'HT4'},{'HTT'},{'target'}];

mydata = table(v2a,v1a,v1b,v4,vT,target,'VariableNames',varnames);

input_filename = 'results/dominance/LZcorr_pcb';
writetable(mydata, [basedir,input_filename,note, '.csv']);

%call from terminal:
% python3 ~/Documents/GIT/DMT_NCT/fxns/dominance_analysis_al857.py ~/Documents/GIT/DMT_NCT/results/dominance/LZcorr_pcb_gsr_volnorm.csv  target 5 0 ~/Documents/GIT/DMT_NCT/results/dominance/LZcorr_pcb_gsr_volnorm


%% regional CE vs intensity

load(fullfile([basedir,'results/regional_continuous_CE_PCB',note,'.mat']),'global_CE_pcb','regional_CE_pcb');

load([basedir,'data/intensity_ratings.mat'])

m_dmt_int = mean(dmt_intensity);
m_pcb_int = mean(pcb_intensity);

load([basedir,'data/5HTvecs_sch116.mat'], 'mean5HT2A_sch116')
receptor_vec = zscore(mean5HT2A_sch116);

TR=2; window=60/TR;

win_dmt_ce=[];
win_pcb_ce=[];

%get avg energy for each min
k=1; 
for i=1:length(m_pcb_int)
    if i==length(m_pcb_int)
        win_pcb_ce(:,i,:) = squeeze(nanmean(regional_CE_pcb(:,k:end,:),2));

    else
        win_pcb_ce(:,i,:) = squeeze(nanmean(regional_CE_pcb(:,k:k+window,:),2));

    end
    k=k+window;
end

%% corr with sub

[r_avg,p_avg]=corr(m_pcb_int',squeeze(mean(win_pcb_ce,1)),'type','Spearman');

[r2a,p2a] = corr(r_avg',receptor_vec,'type','Spearman')
figure;
scatter(receptor_vec,r_avg'); lsline
xlabel('2a density')
ylabel('Rho: Nodal correlation b/w CE and intensity')
text(-2,-0.5,['rho = ',num2str(r2a),'; p = ',num2str(p2a)])

%% FIGURE 3a, right.  gummibrain from Keith Jamison-> https://github.com/kjamison/atlasblobs
load(fullfile([basedir,'fxns/gummi_kj/atlasblobs_saved.mat']))

%choose atlas
whichatlas={'sch116'};

%set colormap
cmap = brewermap([],'*RdBu');

data = r_avg';
clim=[-0.8 0.8];
 
f=figure;

img=display_atlas_blobs(data,atlasblobs,...
    'atlasname',whichatlas,...
    'render',true,...
    'backgroundimage',true,...
    'crop',true,...
    'colormap',cmap,...
    'clim', clim,...
    'backgroundimage','none',...
    'backgroundcolor',[1 1 1]);%,...
%         'roimask',roimask); %last argmument optional, depends if you need to remove roi's

imshow(img);
c=colorbar('SouthOutside', 'fontsize', 16);
c.Label.String='Pearson''s R';
set(gca,'colormap',cmap);
caxis(clim);

%% FIGURE 3c, right (above with no sub and spin test)

load([basedir,'fxns/SpinTests/rotated_maps/rotated_Schaefer_100.mat'])

r_avg=r_avg(1:100);
rec_nosub = receptor_vec(1:100);

%true distribution correlation, will compare strength of against spins
[r2a,p2a] = corr(r_avg',rec_nosub,'type','Spearman')

%permtest
p_spin_vasa_right = perm_sphere_p_al857(r_avg',rec_nosub,perm_id,'Spearman')


figure;
scatter(rec_nosub,r_avg'); lsline
xlabel('2a density')
ylabel('Rho: Nodal correlation b/w CE and intensity')
text(-.7,0.1,['rho = ',num2str(r2a),'; p = ',num2str(p_spin_vasa_right)])

%% dominance analysis - see figure4.m for plotting these data
load([basedir,'data/5HTvecs_sch116.mat'])

v2a=mean5HT2A_sch116;

v1a=mean5HT1A_sch116;

v1b=mean5HT1B_sch116;

v4=mean5HT4_sch116;

vT=mean5HTT_sch116;

allv = [v2a v1a v1b v4 vT];

[r_avg,p_avg]=corr(m_pcb_int',squeeze(mean(win_pcb_ce,1)));
target = r_avg';

varnames=[{'HT2a'},{'HT1a'},{'HT1b'},{'HT4'},{'HTT'},{'target'}];

mydata = table(v2a,v1a,v1b,v4,vT,target,'VariableNames',varnames);

input_filename = 'results/dominance/intcorr_pcb';
writetable(mydata, [basedir,input_filename,note, '.csv']);

%call from terminal:
% python3 ~/Documents/GIT/DMT_NCT/fxns/dominance_analysis_al857.py ~/Documents/GIT/DMT_NCT/results/dominance/intcorr_pcb_gsr_volnorm.csv  target 5 0 ~/Documents/GIT/DMT_NCT/results/dominance/intcorr_pcb_gsr_volnorm




%% correct

pfdr = mafdr([p_spin_vasa_left p_spin_vasa_mid p_spin_vasa_right],'BHFDR', 1)