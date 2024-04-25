%%  Load time-resolved CE from <gen_time_resolved_ce.m> outputs 
%   - Generate Figure 3 and related SI figures
%   - Additionally compile and compute Dominance statistics (Plotted in
%   <figure4.m> script.

clear all; close all;

basedir = '~/Documents/GIT/DMT_NCT/';

note='_gsr_volnorm'; %track different proc streams

nsub=14;

%% SI FIGURE correlation including subcortex, left

load(fullfile([basedir,'results/regional_continuous_prepost_CE_DMT',note,'.mat']),'global_CE_dmt','regional_CE_dmt','global_CE_dmt_pre','regional_CE_dmt_pre');

load([basedir,'data/5HTvecs_sch116.mat'], 'mean5HT2A_sch116');
receptor_vec = zscore(mean5HT2A_sch116);

pre_energy = squeeze(mean(regional_CE_dmt_pre,2));
post_energy = squeeze(mean(regional_CE_dmt,2));

diff_energy = (post_energy-pre_energy)./(pre_energy);

[r2a_sub_left,p2a_sub_left] = corr(mean(diff_energy)',receptor_vec,'type','Spearman')


figure;
scatter(receptor_vec,mean(diff_energy)'); lsline
xlabel('2a density')
ylabel('pre/post DMT Nodal Control Energy (∆)')
text(-2,-0.4,['rho = ',num2str(r2a_sub_left),'; p = ',num2str(p2a_sub_left)])

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
ylabel('pre/post DMT Nodal Control Energy (∆)')
text(-.7,-0.5,['rho = ',num2str(r2a),'; p = ',num2str(p_spin_vasa_left)])


%% dominance analysis - see figure4.m for plotting these data
load([basedir,'data/5HTvecs_sch116.mat'])

pre_energy = squeeze(mean(regional_CE_dmt_pre,2));
post_energy = squeeze(mean(regional_CE_dmt,2));

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

input_filename = 'results/dominance/diffE_continuous';
writetable(mydata, [basedir,input_filename,note, '.csv']);

%call from terminal:
% python3 ~/Documents/GIT/DMT_NCT/fxns/dominance_analysis_al857.py ~/Documents/GIT/DMT_NCT/results/dominance/diffE_continuous_gsr_volnorm.csv  target 5 0 ~/Documents/GIT/DMT_NCT/results/dominance/diffE_continuous_gsr_volnorm




%% SI FIGURE correlation including subcortex, middle, correlate nodal CE with LZ

load(fullfile([basedir,'results/regional_continuous_CE_DMT',note,'.mat']),'global_CE_dmt','regional_CE_dmt');

load([basedir,'data/RegressorLZInterpscrubbedConvolvedAvg.mat'])

m_dmt_LZ = nanmean(RegDMT2);
m_dmt_LZ = m_dmt_LZ(2:839); %first frame is excluded from CE calcs


[r_avg,p_avg]=corr(m_dmt_LZ',squeeze(mean(regional_CE_dmt,1)),'type','Spearman');


[r2a_sub_middle,p2a_sub_middle] = corr(r_avg',receptor_vec,'type','Spearman')
figure;
scatter(receptor_vec,r_avg'); lsline
xlabel('2a density')
ylabel('Rho: Nodal correlation b/w CE and LZ')
text(-2,-0.5,['rho = ',num2str(r2a_sub_middle),'; p = ',num2str(p2a_sub_middle)])

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

p_spin_vasa_middle = perm_sphere_p_al857(r_avg',rec_nosub,perm_id,'Spearman')

figure;
scatter(rec_nosub,r_avg'); lsline
xlabel('2a density')
ylabel('Rho: Nodal correlation b/w CE and LZ')
text(-.7,-0.5,['rho = ',num2str(r2a),'; p = ',num2str(p_spin_vasa_middle)])

%% dominance analysis - see figure4.m for plotting these data
load([basedir,'data/5HTvecs_sch116.mat'])

v2a=mean5HT2A_sch116;

v1a=mean5HT1A_sch116;

v1b=mean5HT1B_sch116;

v4=mean5HT4_sch116;

vT=mean5HTT_sch116;

allv = [v2a v1a v1b v4 vT];

[r_avg,p_avg]=corr(m_dmt_LZ',squeeze(mean(regional_CE_dmt,1)));
target = r_avg';

varnames=[{'HT2a'},{'HT1a'},{'HT1b'},{'HT4'},{'HTT'},{'target'}];

mydata = table(v2a,v1a,v1b,v4,vT,target,'VariableNames',varnames);

input_filename = 'results/dominance/LZcorr';
writetable(mydata, [basedir,input_filename,note, '.csv']);

%call from terminal:
% python3 ~/Documents/GIT/DMT_NCT/fxns/dominance_analysis_al857.py ~/Documents/GIT/DMT_NCT/results/dominance/LZcorr_gsr_volnorm.csv  target 5 0 ~/Documents/GIT/DMT_NCT/results/dominance/LZcorr_gsr_volnorm

%% regional CE vs intensity

load(fullfile([basedir,'results/regional_continuous_CE_DMT',note,'.mat']),'global_CE_dmt','regional_CE_dmt');

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
for i=1:length(m_dmt_int)
    if i==length(m_dmt_int)
        win_dmt_ce(:,i,:) = squeeze(nanmean(regional_CE_dmt(:,k:end,:),2));

    else
        win_dmt_ce(:,i,:) = squeeze(nanmean(regional_CE_dmt(:,k:k+window,:),2));

    end
    k=k+window;
end

%% SI FIGURE correlation including subcortex, right

[r_avg,p_avg]=corr(m_dmt_int',squeeze(mean(win_dmt_ce,1)),'type','Spearman');

[r2a_sub_right,p2a_sub_right] = corr(r_avg',receptor_vec,'type','Spearman')
figure;
scatter(receptor_vec,r_avg'); lsline
xlabel('2a density')
ylabel('Rho: Nodal correlation b/w CE and intensity')
text(-2,-0.5,['rho = ',num2str(r2a_sub_right),'; p = ',num2str(p2a_sub_right)])

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
text(-.7,-0.8,['rho = ',num2str(r2a),'; p = ',num2str(p_spin_vasa_right)])

%% dominance analysis - see figure4.m for plotting these data
load([basedir,'data/5HTvecs_sch116.mat'])

v2a=mean5HT2A_sch116;

v1a=mean5HT1A_sch116;

v1b=mean5HT1B_sch116;

v4=mean5HT4_sch116;

vT=mean5HTT_sch116;

allv = [v2a v1a v1b v4 vT];

[r_avg,p_avg]=corr(m_dmt_int',squeeze(mean(win_dmt_ce,1)));
target = r_avg';

varnames=[{'HT2a'},{'HT1a'},{'HT1b'},{'HT4'},{'HTT'},{'target'}];

mydata = table(v2a,v1a,v1b,v4,vT,target,'VariableNames',varnames);

input_filename = 'results/dominance/intcorr';
writetable(mydata, [basedir,input_filename,note, '.csv']);

%call from terminal:
% python3 ~/Documents/GIT/DMT_NCT/fxns/dominance_analysis_al857.py ~/Documents/GIT/DMT_NCT/results/dominance/intcorr_gsr_volnorm.csv  target 5 0 ~/Documents/GIT/DMT_NCT/results/dominance/intcorr_gsr_volnorm

%% correct pvalues

pfdr_main = mafdr([p_spin_vasa_left p_spin_vasa_middle p_spin_vasa_right],'BH',1)
pfdr_si = mafdr([p2a_sub_left p2a_sub_middle p2a_sub_right],'BH',1)