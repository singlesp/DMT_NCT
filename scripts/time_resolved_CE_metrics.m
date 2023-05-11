%% tr by tr E of the continuous DMT scans - (FIGURE 3, middle, right)
clear all; close all;

basedir = '~/Documents/GIT/DMT_NCT/';
load([basedir,'data/DMT_clean_mni_continuous_fullPreprocsch116.mat'],'ts_gsr')
load([basedir,'data/Schaefer116_HCP_DTI_count.mat'],'vol_normalized_sc')


note='_gsr_volnorm'; %track different proc streams
sc=vol_normalized_sc; 
TS = ts_gsr; %change depending on proc stream

nsub=14;

[nparc, transitions] = size(TS{1,1});
transitions = transitions - 2; %skip first transition

%% regional CE DMT

global_CE_dmt = NaN(nsub,transitions);
regional_CE_dmt = NaN(nsub,transitions,nparc);

T=1; %time-horizon
c=1; %adjacency matrix normalization factor

InputVector = ones(nparc,1);
B = InputVector .*eye(nparc) + eye(nparc); 
Anorm = NORMALIZE(sc,c);

for i = 1:nsub

    ts = TS{i,1};
    ts(:,1) = []; %discard first frame
    
    [global_CE_dmt(i,:),regional_CE_dmt(i,:,:)] = time_resolved_control_energy(Anorm,T,B,ts);

end

save(fullfile([basedir,'results/regional_continuous_CE_DMT',note,'.mat']),'global_CE_dmt','regional_CE_dmt');


%% regional CE PCB

global_CE_pcb = NaN(nsub,transitions);
regional_CE_pcb = NaN(nsub,transitions,nparc);

for i = 1:nsub

    ts = TS{i,2};
    ts(:,1) = []; %discard first frame
    
    [global_CE_pcb(i,:),regional_CE_pcb(i,:,:)] = time_resolved_control_energy(Anorm,T,B,ts);

end

save(fullfile([basedir,'results/regional_continuous_CE_PCB',note,'.mat']),'global_CE_pcb','regional_CE_pcb');



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

[r_avg,p_avg]=corr(m_dmt_int',squeeze(mean(win_dmt_ce,1)));

%% SI Figure 4, right
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
p_spin_vasa = perm_sphere_p_al857(r_avg',rec_nosub,perm_id,'Spearman')


figure;
scatter(rec_nosub,r_avg'); lsline
xlabel('2a density')
ylabel('Rho: Nodal correlation b/w CE and intensity')
text(-.7,-0.8,['rho = ',num2str(r2a),'; p = ',num2str(p_spin_vasa)])

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


%% SI FIGURE 4, middle correlate nodal CE with LZ

load([basedir,'data/RegressorLZInterpscrubbedConvolvedAvg.mat'])

m_dmt_LZ = nanmean(RegDMT2);
m_dmt_LZ = m_dmt_LZ(2:839); %first frame is excluded from CE calcs


[r_avg,p_avg]=corr(m_dmt_LZ',squeeze(mean(regional_CE_dmt,1)));


[r2a,p2a] = corr(r_avg',receptor_vec,'type','Spearman')
figure;
scatter(receptor_vec,r_avg'); lsline
xlabel('2a density')
ylabel('Rho: Nodal correlation b/w CE and LZ')
text(-2,-0.5,['rho = ',num2str(r2a),'; p = ',num2str(p2a)])

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

p_spin_vasa = perm_sphere_p_al857(r_avg',rec_nosub,perm_id,'Spearman')

figure;
scatter(rec_nosub,r_avg'); lsline
xlabel('2a density')
ylabel('Rho: Nodal correlation b/w CE and LZ')
text(-.7,-0.5,['rho = ',num2str(r2a),'; p = ',num2str(p_spin_vasa)])

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
