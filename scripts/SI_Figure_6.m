%% SI FIGURE 6

clear all; close all;

basedir = '~/Documents/GIT/DMT_NCT/';
load([basedir,'data/DMT_clean_mni_continuous_fullPreprocsch116.mat'],'ts_gsr')
load([basedir,'data/Schaefer116_HCP_DTI_count.mat'],'vol_normalized_sc')


note='_gsr_volnorm'; %track different proc streams
sc=vol_normalized_sc; 
TS = ts_gsr; %change depending on proc stream

nsub=14;

T=1; c=1;
Anorm = NORMALIZE(sc,c); 
WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T

%% uses a faster version of control energy calculation - can be used when B matrix is uniform - values will be equivalent to the slow version, albeit scaled


for i = 1:nsub

    ts = TS{i,1}; %full continuous dmt scan

    %% no normalization
    ts_pre = ts(:,1:240); ts_post = ts(:,241:480);
    ts_pre(:,1) = []; ts_post(:,1) = []; %discard first frame
    
    
    E_pre_dmt(i,:) = time_resolved_control_energy_fast(Anorm,T,ts_pre);
    E_tot_pre_dmt(i,1) = mean(E_pre_dmt(i,:),2);
   
    E_post_dmt(i,:) = time_resolved_control_energy_fast(Anorm,T,ts_post);
    E_tot_post_dmt(i,1) = mean(E_post_dmt(i,:),2);
    
    %% L2 normalization
    
    ts_pre = L2MAGNITUDENORM(ts(:,1:240)); ts_pre(:,1)=[];
    ts_post = L2MAGNITUDENORM(ts(:,241:480)); ts_post(:,1)=[];
    
    E_pre_dmt_l2(i,:) = time_resolved_control_energy_fast(Anorm,T,ts_pre);
    E_tot_pre_dmt_l2(i,1) = mean(E_pre_dmt_l2(i,:),2);
   
    E_post_dmt_l2(i,:) = time_resolved_control_energy_fast(Anorm,T,ts_post);
    E_tot_post_dmt_l2(i,1) = mean(E_post_dmt_l2(i,:),2);

    %% Distance normalization - requires some special handling - must break into pairs of states prior to normalization
    ts_pre = ts(:,1:240); ts_post = ts(:,241:480);
    
    x0_pre = ts_pre(:,2:size(ts_pre,2)-1); %% discard first volume
    xf_pre = ts_pre(:,3:size(ts_pre,2));
    
    [x0_pre,xf_pre,x0xfmag_pre(i,:)] = DISTANCENORM(x0_pre,xf_pre);
    
    x0_post = ts_post(:,2:size(ts_post,2)-1); %% discard first volume
    xf_post = ts_post(:,3:size(ts_post,2));
    
    [x0_post,xf_post,x0xfmag_post(i,:)] = DISTANCENORM(x0_post,xf_post);
    
    E_pre_dmt_dist(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0_pre, xf_pre, T,false);
    E_tot_pre_dmt_dist(i,1) = mean(E_pre_dmt_dist(i,:),2);
   
    E_post_dmt_dist(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0_post, xf_post, T,false);
    E_tot_post_dmt_dist(i,1) = mean(E_post_dmt_dist(i,:),2);
    
    %% Double normalization - requires some special handling - must break into pairs of states prior to normalization
    x0_pre = ts_pre(:,2:size(ts_pre,2)-1); %% discard first volume
    xf_pre = ts_pre(:,3:size(ts_pre,2));
    
    [x0_pre,xf_pre,x0xfmag_pre_double(i,:)] = DOUBLENORM(x0_pre,xf_pre);
    
    x0_post = ts_post(:,2:size(ts_post,2)-1); %% discard first volume
    xf_post = ts_post(:,3:size(ts_post,2));
    
    [x0_post,xf_post,x0xfmag_post_double(i,:)] = DOUBLENORM(x0_post,xf_post);
    
    E_pre_dmt_double(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0_pre, xf_pre, T,false);
    E_tot_pre_dmt_double(i,1) = mean(E_pre_dmt_double(i,:),2);
   
    E_post_dmt_double(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0_post, xf_post, T,false);
    E_tot_post_dmt_double(i,1) = mean(E_post_dmt_double(i,:),2);
    
    %% radial normalization - requires some special handling - must break into pairs of states prior to normalization
    x0_pre = ts_pre(:,2:size(ts_pre,2)-1); %% discard first volume
    xf_pre = ts_pre(:,3:size(ts_pre,2));
    
    [x0_pre,xf_pre] = RADIALNORM(x0_pre,xf_pre);
    
    x0_post = ts_post(:,2:size(ts_post,2)-1); %% discard first volume
    xf_post = ts_post(:,3:size(ts_post,2));
    
    [x0_post,xf_post] = RADIALNORM(x0_post,xf_post);
    
    E_pre_dmt_rad(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0_pre, xf_pre, T,false);
    E_tot_pre_dmt_rad(i,1) = mean(E_pre_dmt_rad(i,:),2);
   
    E_post_dmt_rad(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0_post, xf_post, T,false);
    E_tot_post_dmt_rad(i,1) = mean(E_post_dmt_rad(i,:),2);
    
end


%%
norm_types = [{'None'}; {'Radial'}; {'L2'}; {'Distance'}; {'Double'}];
pre_dmt = [E_tot_pre_dmt E_tot_pre_dmt_rad E_tot_pre_dmt_l2 E_tot_pre_dmt_dist E_tot_pre_dmt_double];
post_dmt = [E_tot_post_dmt E_tot_post_dmt_rad E_tot_post_dmt_l2 E_tot_post_dmt_dist E_tot_post_dmt_double];


[h,p,~,t] = ttest(pre_dmt,post_dmt);
pfdr = mafdr(p,'BH',1);


figure;
for i=1:5
    if i==1
        subplot(3,2,1.5)
    else
        subplot(3,2,i+1)
    end
    data = [pre_dmt(:,i) post_dmt(:,i)];
    violin(data)
    text(1.5,0.95*max(max(data)),[[{'t = '}, {num2str(t.tstat(i))}];[{'p = '},{num2str(pfdr(i))}]]);
    title(char(norm_types{i}))
    xticks([1 2])
    xticklabels([{'preDMT'},{'postDMT'}])
end
%%
figure;
subplot(1,2,1)
data = [mean(x0xfmag_pre,2) mean(x0xfmag_post,2)];
violin(data)
[~,p,~,t] = ttest(data(:,1),data(:,2))
title('inter-state distance between un-normalized states')
text(1.5,0.95*max(max(data)),[[{'t = '}, {num2str(t.tstat)}];[{'p = '},{num2str(pfdr(i))}]]);
xticks([1 2])
xticklabels([{'preDMT'},{'postDMT'}])

subplot(1,2,2)
data = [mean(x0xfmag_pre_double,2) mean(x0xfmag_post_double,2)];
violin(data)
[~,p,~,t] = ttest(data(:,1),data(:,2))
title('inter-state distance between L2-normalized states')
text(1.5,0.95*max(max(data)),[[{'t = '}, {num2str(t.tstat)}];[{'p = '},{num2str(pfdr(i))}]]);
xticks([1 2])
xticklabels([{'preDMT'},{'postDMT'}])

%%
y = mean(x0xfmag_post,2)-mean(x0xfmag_pre,2);
x = mean(x0xfmag_post_double,2)-mean(x0xfmag_pre_double,2);
scatter_corr(x,y,[{'avg dist b/w L2-normalized states (post-pre)'},{'avg dist b/w un-normalized states (post-pre)'}])


y = E_tot_post_dmt-E_tot_pre_dmt;
x = E_tot_post_dmt_l2-E_tot_pre_dmt_l2;
scatter_corr(x,y,[{'avg energy w/ L2-normalized states (post-pre)'},{'avg energy w/ un-normalized states (post-pre)'}])

%% 
y = mean(x0xfmag_post,2)-mean(x0xfmag_pre,2);
x = E_tot_post_dmt-E_tot_pre_dmt;
scatter_corr(x,y,[{'avg energy w/ un-normalized states (post-pre)'},{'avg dist b/w un-normalized states (post-pre)'}])

y = mean(x0xfmag_post_double,2)-mean(x0xfmag_pre_double,2);
x = E_tot_post_dmt_l2-E_tot_pre_dmt_l2;
scatter_corr(x,y,[{'avg energy w/ L2-normalized states (post-pre)'},{'avg dist b/w L2-normalized states (post-pre)'}])
