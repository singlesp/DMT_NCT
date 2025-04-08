%% Get time-resolved control energy. First full scans, then pre/post 8 minutes
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

%% now pre/post

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
