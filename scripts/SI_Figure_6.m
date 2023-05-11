%% SI FIGURE 6

clear all; close all;

basedir = '~/Documents/GIT/DMT_NCT/';
load([basedir,'data/DMT_clean_mins8_mni_sch116.mat'],'ts_gsr')
load([basedir,'data/Schaefer116_HCP_DTI_count.mat'], 'vol_normalized_sc')

note='_gsr_volnorm'; %track different proc streams
sc=vol_normalized_sc;
TS = ts_gsr; %change depending on proc stream

nsub=14;

T=1; c=1;
Anorm = NORMALIZE(sc,c); 
WcI = GRAMIAN_FAST(Anorm, T); % compute gramian inverse for control horizon T

%% uses a faster version of control energy calculation - can be used when B matrix is uniform - values will be equivalent to the slow version, albeit scaled


for i = 1:nsub

    %% no normalization
    ts_pre = TS{i,1};
    ts_post = TS{i,2};
    
    x0_pre = ts_pre(:,2:size(ts_pre,2)-1); %% discard first volume
    xf_pre = ts_pre(:,3:size(ts_pre,2));
    
    x0_post = ts_post(:,2:size(ts_post,2)-1);
    xf_post = ts_post(:,3:size(ts_post,2));
    
    E_pre_dmt(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0_pre, xf_pre, T,false);
    E_tot_pre_dmt(i,1) = mean(E_pre_dmt(i,:)');
   
    E_post_dmt(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0_post, xf_post, T,false);
    E_tot_post_dmt(i,1) = mean(E_post_dmt(i,:)');
    
    %% L2 normalization
    
    ts_pre = L2MAGNITUDENORM(TS{i,1});
    ts_post = L2MAGNITUDENORM(TS{i,2});
    
    x0_pre = ts_pre(:,2:size(ts_pre,2)-1); %% discard first volume
    xf_pre = ts_pre(:,3:size(ts_pre,2));
    
    x0_post = ts_post(:,2:size(ts_post,1)-1);
    xf_post = ts_post(:,3:size(ts_post,1));
    
    E_pre_dmt_l2(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0_pre, xf_pre, T,false);
    E_tot_pre_dmt_l2(i,1) = mean(E_pre_dmt_l2(i,:)');
   
    E_post_dmt_l2(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0_post, xf_post, T,false);
    E_tot_post_dmt_l2(i,1) = mean(E_post_dmt_l2(i,:)');

    %% Distance normalization
    ts_pre = TS{i,1};
    ts_post = TS{i,2};
    
    x0_pre = ts_pre(:,2:size(ts_pre,2)-1); %% discard first volume
    xf_pre = ts_pre(:,3:size(ts_pre,2));
    
    [x0_pre,xf_pre,x0xfmag_pre(i,:)] = DISTANCENORM(x0_pre,xf_pre);
    
    x0_post = ts_post(:,2:size(ts_post,2)-1);
    xf_post = ts_post(:,3:size(ts_post,2));
    
    [x0_post,xf_post,x0xfmag_post(i,:)] = DISTANCENORM(x0_post,xf_post);
    
    E_pre_dmt_dist(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0_pre, xf_pre, T,false);
    E_tot_pre_dmt_dist(i,1) = mean(E_pre_dmt_dist(i,:)');
   
    E_post_dmt_dist(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0_post, xf_post, T,false);
    E_tot_post_dmt_dist(i,1) = mean(E_post_dmt_dist(i,:)');
    
    %% Double normalization
    x0_pre = ts_pre(:,2:size(ts_pre,2)-1); %% discard first volume
    xf_pre = ts_pre(:,3:size(ts_pre,2));
    
    [x0_pre,xf_pre,x0xfmag_pre_double(i,:)] = DOUBLENORM(x0_pre,xf_pre);
    
    x0_post = ts_post(:,2:size(ts_post,2)-1);
    xf_post = ts_post(:,3:size(ts_post,2));
    
    [x0_post,xf_post,x0xfmag_post_double(i,:)] = DOUBLENORM(x0_post,xf_post);
    
    E_pre_dmt_double(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0_pre, xf_pre, T,false);
    E_tot_pre_dmt_double(i,1) = mean(E_pre_dmt_double(i,:)');
   
    E_post_dmt_double(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0_post, xf_post, T,false);
    E_tot_post_dmt_double(i,1) = mean(E_post_dmt_double(i,:)');
    
    %% radial normalization
    x0_pre = ts_pre(:,2:size(ts_pre,2)-1); %% discard first volume
    xf_pre = ts_pre(:,3:size(ts_pre,2));
    
    [x0_pre,xf_pre] = RADIALNORM(x0_pre,xf_pre);
    
    x0_post = ts_post(:,2:size(ts_post,2)-1);
    xf_post = ts_post(:,3:size(ts_post,2));
    
    [x0_post,xf_post] = RADIALNORM(x0_post,xf_post);
    
    E_pre_dmt_rad(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0_pre, xf_pre, T,false);
    E_tot_pre_dmt_rad(i,1) = mean(E_pre_dmt_rad(i,:)');
   
    E_post_dmt_rad(i,:) = MIN_CONTROL_ENERGY(Anorm, WcI, x0_post, xf_post, T,false);
    E_tot_post_dmt_rad(i,1) = mean(E_post_dmt_rad(i,:)');
    
end


%%
norm_types = [{'None'}; {'Radial'}; {'L2'}; {'Distance'}; {'Double'}];
pre_dmt = [E_tot_pre_dmt E_tot_pre_dmt_rad E_tot_pre_dmt_l2 E_tot_pre_dmt_dist E_tot_pre_dmt_double];
post_dmt = [E_tot_post_dmt E_tot_post_dmt_rad E_tot_post_dmt_l2 E_tot_post_dmt_dist E_tot_post_dmt_double];


[h,p,~,t] = ttest(pre_dmt,post_dmt);


figure;
for i=1:5
    if i==1
        subplot(3,2,1.5)
    else
        subplot(3,2,i+1)
    end
    data = [pre_dmt(:,i) post_dmt(:,i)];
    violin(data)
    text(1.5,0.95*max(max(data)),[[{'t = '}, {num2str(t.tstat(i))}];[{'p = '},{num2str(p(i))}]]);
    title(char(norm_types{i}))
    xticks([1 2])
    xticklabels([{'preDMT'},{'postDMT'}])
end

figure;
subplot(1,2,1)
violin([mean(x0xfmag_pre,1)' mean(x0xfmag_post,1)'])
title('inter-state distance between un-normalized states')
xticks([1 2])
xticklabels([{'preDMT'},{'postDMT'}])

subplot(1,2,2)
violin([mean(x0xfmag_pre_double,1)' mean(x0xfmag_post_double,1)'])
title('inter-state distance between L2-normalized states')
xticks([1 2])
xticklabels([{'preDMT'},{'postDMT'}])

