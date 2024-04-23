%% regional CE for mins8 (Figure 3, left)
clear all; close all;

basedir = '~/Documents/GIT/DMT_NCT/';
load([basedir,'data/DMT_clean_mni_continuous_fullPreprocsch116.mat'],'ts_gsr')
load([basedir,'data/Schaefer116_HCP_DTI_count.mat'], 'vol_normalized_sc')

sc=vol_normalized_sc;

note='_gsr_volnorm'; %track different proc streams
TS = ts_gsr; %change depending on proc stream

nsub=14;


[nparc, ~] = size(TS{1,1});

transitions = 240 - 2; %first tr skipped

%% regional CE DMT

global_CE_dmt = NaN(nsub,transitions);
regional_CE_dmt = NaN(nsub,transitions,nparc);
global_CE_dmt_pre = NaN(nsub,transitions);
regional_CE_dmt_pre = NaN(nsub,transitions,nparc);

T=1; %time-horizon
c=1; %adjacency matrix normalization factor

InputVector = ones(nparc,1);
B = InputVector .*eye(nparc) + eye(nparc); 
Anorm = NORMALIZE(sc,c);

for i = 1:nsub

    ts = TS{i,1};
    ts_pre = ts(:,1:240); ts_post = ts(:,241:480);
    ts_post(:,1) = []; %discard first frame
    
    [global_CE_dmt(i,:),regional_CE_dmt(i,:,:)] = time_resolved_control_energy(Anorm,T,B,ts_post);
    
    ts_pre(:,1) = []; %discard first frame
    
    [global_CE_dmt_pre(i,:),regional_CE_dmt_pre(i,:,:)] = time_resolved_control_energy(Anorm,T,B,ts_pre);

end

save(fullfile([basedir,'results/regional_continuous_prepost_CE_DMT',note,'.mat']),'global_CE_dmt','regional_CE_dmt','global_CE_dmt_pre','regional_CE_dmt_pre');


%% regional CE PCB


global_CE_pcb = NaN(nsub,transitions);
regional_CE_pcb = NaN(nsub,transitions,nparc);
global_CE_pcb_pre = NaN(nsub,transitions);
regional_CE_pcb_pre = NaN(nsub,transitions,nparc);

T=1; %time-horizon
c=1; %adjacency matrix normalization factor

InputVector = ones(nparc,1);
B = InputVector .*eye(nparc) + eye(nparc); 
Anorm = NORMALIZE(sc,c);

for i = 1:nsub

    ts = TS{i,2};
    ts_pre = ts(:,1:240); ts_post = ts(:,241:480);
    ts_post(:,1) = []; %discard first frame
    
    [global_CE_pcb(i,:),regional_CE_pcb(i,:,:)] = time_resolved_control_energy(Anorm,T,B,ts_post);
    
    ts_pre(:,1) = []; %discard first frame
    
    [global_CE_pcb_pre(i,:),regional_CE_pcb_pre(i,:,:)] = time_resolved_control_energy(Anorm,T,B,ts_pre);


end

save(fullfile([basedir,'results/regional_continuous_prepost_CE_PCB',note,'.mat']),'global_CE_pcb','regional_CE_pcb','global_CE_pcb_pre','regional_CE_pcb_pre');


%% SI FIGURE 4, left

load(fullfile([basedir,'results/regional_continuous_prepost_CE_DMT',note,'.mat']),'global_CE_dmt','regional_CE_dmt','global_CE_dmt_pre','regional_CE_dmt_pre');

load([basedir,'data/5HTvecs_sch116.mat'], 'mean5HT2A_sch116');
receptor_vec = zscore(mean5HT2A_sch116);

pre_energy = squeeze(mean(regional_CE_dmt_pre,2));
post_energy = squeeze(mean(regional_CE_dmt,2));

diff_energy = (post_energy-pre_energy)./(pre_energy);

[r2a,p2a] = corr(mean(diff_energy)',receptor_vec,'type','Spearman')


figure;
scatter(receptor_vec,mean(diff_energy)'); lsline
xlabel('2a density')
ylabel('pre/post DMT Nodal Control Energy (∆)')
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

p_spin_vasa = perm_sphere_p_al857(mean(diff_energy)',rec_nosub,perm_id,'Spearman')

figure;
scatter(rec_nosub,mean(diff_energy)'); lsline
xlabel('2a density')
ylabel('pre/post DMT Nodal Control Energy (∆)')
text(-.7,-0.5,['rho = ',num2str(r2a),'; p = ',num2str(p_spin_vasa)])


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
