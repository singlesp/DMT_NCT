%% main script for simulated E

clear all; close all;

basedir = '~/Documents/GIT/DMT_NCT/';

note='_gsr_volnorm'; % track separate proc streams

load([basedir,'results/regional_continuous_CE_DMT',note,'.mat'])
load([basedir,'results/regional_continuous_CE_PCB',note,'.mat'])

load([basedir,'data/DMT_clean_mni_continuous_fullPreprocsch116.mat'],'ts_gsr')
load([basedir,'data/Schaefer116_HCP_DTI_count.mat'], 'vol_normalized_sc')

TS = ts_gsr; %change depending on proc streams
sc=vol_normalized_sc;

nsub=14;

[nparc, transitions] = size(TS{1,1});
transitions = transitions - 2; %skip first transition

addpath(genpath('~/Documents/MATLAB/boundedline'))

%% simulated E using Alpha effect compartment concentration

load data/5HTvecs_sch116.mat mean5HT2A_sch116
scaled2a = mean5HT2A_sch116./max(mean5HT2A_sch116);

sim_conc =  readtable([basedir,'/data/simulated_effectcomp_concentrations.csv']);
dmt_conc = [zeros(239,1); double(sim_conc.CE_alpha(17:615))];

T=1;
c=1;
Anorm = NORMALIZE(sc,c);

E_weighted_regional_alpha = NaN(nsub,transitions,nparc);
E_weighted_global_alpha = NaN(nsub,transitions);

for rho=[1 10 20 30 40 50 60 70]
    
    for i=1:transitions
        InputVector(:,i) = 1 + rho*(dmt_conc(i)*scaled2a);
    end
    
    for i = 1:nsub
        
        ts = TS{i,2};

        
        x0 = ts(:,2:size(ts,2)-1);
        xf = ts(:,3:size(ts,2));

        for transition = 1:transitions
            B = InputVector(:,transition) .*eye(nparc) + eye(nparc);
            [x, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
            E_weighted_regional_alpha(i,transition,:) = sum(u.^2)*T/1001;
            E_weighted_global_alpha(i,transition) = sum(sum(u.^2))*T/1001; % integrate over inputs
        end
        
    end
    
    save(fullfile([basedir,'results/simulations/simulated_sub_Alphaeffectconc_E_rho_',num2str(rho),note,'.mat']),'E_weighted_global_alpha','E_weighted_regional_alpha');

end

%%

mean_dmt = nanmean(global_CE_dmt);
mean_pcb = nanmean(global_CE_pcb);

i=1;
for rho=[1 10 20 30 40 50 60 70]
    
    load(fullfile([basedir,'results/simulations/simulated_sub_Alphaeffectconc_E_rho_',num2str(rho),note,'.mat']))
    
    mean_alpha_sim = nanmean(E_weighted_global_alpha);

    %euclidean distance
    r_dmt(i) = sqrt(sum((mean_dmt-mean_alpha_sim).^2));
    r_pcb(i) = sqrt(sum((mean_pcb-mean_alpha_sim).^2));

    i=i+1;
end

r_dmt
r_pcb

figure;
hold on
    plot(r_dmt)
    plot(r_pcb)
xlim([0.75 8.25])
ylim([0 2500])
xticks(1:8)
xticklabels([{'1'},{'10'},{'20'},{'30'},{'40'},{'50'},{'60'},{'70'}]);
xlabel('scaling parameter')
ylabel('Euclidean Distance')
title('Euclidean distance: AlphaEffectConc')
legend('DMT vs Simulation','PCB vs Simulation')

%%

rho=30;
load(fullfile([basedir,'results/simulations/simulated_sub_Alphaeffectconc_E_rho_',num2str(rho),note,'.mat']))

figure;
title([{'Control Energy Time-Series'};{'continuous data (group means)'}])
hold on
    plot_bounded_line(global_CE_dmt,[0 0 0])
    plot_bounded_line(E_weighted_global_alpha,[1 0 0])
ylabel('Control Energy')
ylim([0 300])
tics = linspace(0,28,15);
tics = tics*30;
tics(end)=838;
xticks(tics)
xticklabels([{'-8'},{'-6'},{'-4'},{'-2'},{'0'},{'2'},{'4'},{'6'},{'8'},{'10'},{'12'},{'14'},{'16'},{'18'},{'20'}]);
xlabel('Minutes')
legend('Empirical DMT','Simulated DMT')

%% simulated E using plasma concentration

load data/5HTvecs_sch116.mat mean5HT2A_sch116

scaled2a = mean5HT2A_sch116./max(mean5HT2A_sch116);

sim_conc =  readtable([basedir,'/data/simulated_effectcomp_concentrations.csv']);
dmt_conc = [zeros(239,1); double(sim_conc.plasmaConc(17:615))];

T=1;
c=1;
Anorm = NORMALIZE(sc,c);

for rho=[1 10 20 30 40 50 60 70]
    
    for i=1:transitions
        InputVector(:,i) = 1 + rho*(dmt_conc(i)*scaled2a);
    end
    
    for i = 1:nsub
        
        ts = TS{i,2};
        
        x0 = ts(:,2:size(ts,2)-1);
        xf = ts(:,3:size(ts,2));

        for transition = 1:transitions
            B = InputVector(:,transition) .*eye(nparc) + eye(nparc);
            [x, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
            E_weighted_regional_plasma(i,transition,:) = sum(u.^2)*T/1001;
            E_weighted_global_plasma(i,transition) = sum(sum(u.^2))*T/1001; % integrate over inputs
        end
        
        
    end
    
    save(fullfile([basedir,'results/simulations/simulated_sub_plasmaconc_E_rho_',num2str(rho),note,'.mat']),'E_weighted_global_plasma','E_weighted_regional_plasma');

end

%%

mean_dmt = nanmean(global_CE_dmt);
mean_pcb = nanmean(global_CE_pcb);

i=1;
for rho=[1 10 20 30 40 50 60 70]
    
    load(fullfile([basedir,'results/simulations/simulated_sub_plasmaconc_E_rho_',num2str(rho),note,'.mat']))
    
    mean_plasma_sim = nanmean(E_weighted_global_plasma);
 
    r_dmt(i) = sqrt(sum((mean_dmt-mean_plasma_sim).^2));
    r_pcb(i) = sqrt(sum((mean_pcb-mean_plasma_sim).^2));

    i=i+1;
end

r_dmt
r_pcb

figure;
hold on
    plot(r_dmt)
    plot(r_pcb)
xlim([0.75 8.25])
ylim([0 2500])
xticks(1:8)
xticklabels([{'1'},{'10'},{'20'},{'30'},{'40'},{'50'},{'60'},{'70'}]);
xlabel('scaling parameter')
ylabel('Euclidean Distance')
title('Euclidean distance: PlamaConc')
legend('DMT vs Simulation','PCB vs Simulation')


%% plots for explainer figure; 

% nanmat = NaN(nparc,nparc);
% ind = find(B>0)
% nanmat(ind) = InputVector(:,transition);
% figure; heatmap(nanmat);
% figure; heatmap(nanmat); colormap(iris)
% figure; heatmap(nanmat); colormap(viridis)
% transition=2;
% nanmat(ind) = InputVector(:,transition);
% figure; heatmap(nanmat); colormap(parula)
% figure; heatmap(nanmat); colormap(parula); caxis([1.1 2])

% need to use heatmap and then edit figure properties, change missing data
% color and remove grids

sim_conc =  readtable([basedir,'/data/simulated_effectcomp_concentrations.csv']);
dmt_conc = [zeros(239,1); double(sim_conc.CE_alpha(17:615))];
rho=30;
for i=1:transitions
InputVector(:,i) = 1 + rho*(dmt_conc(i)*scaled2a);
end

figure; plot(InputVector'); ylim([0 3])
figure; plot(InputVector'); ylim([0.5 3])
figure; plot(dmt_conc'); ylim([-0.01 0.07])
figure; imagesc(InputVector)
figure; imagesc(InputVector(:,445)); caxis([min(min(InputVector)) max(max(InputVector))])

%% simulated E using Alpha effect compartment concentration w/ uniform vec

load data/5HTvecs_sch116.mat mean5HT2A_sch116

scaled2a = mean5HT2A_sch116./max(mean5HT2A_sch116);

sim_conc =  readtable([basedir,'data/simulated_effectcomp_concentrations.csv']);
dmt_conc = [zeros(239,1); double(sim_conc.CE_alpha(17:615))];


%this keeps the same amount of total control as true2a, 
% but evenly distributes it across all regions
univec = repelem(sum(scaled2a)/nparc,nparc,1); 

T=1;
c=1;
Anorm = NORMALIZE(sc,c);

for rho=[1 10 20 30 40 50 60 70]
    
    for i=1:transitions
        InputVector(:,i) = 1 + rho*(dmt_conc(i)*univec);
    end
    
    for i = 1:nsub
        
        ts = TS{i,2};
        
        x0 = ts(:,2:size(ts,2)-1);
        xf = ts(:,3:size(ts,2));
       
        for transition = 1:transitions
            B = InputVector(:,transition) .*eye(nparc) + eye(nparc);
            [x, u] = MIN_ENG_CONT(Anorm, T, B, x0(:,transition), xf(:,transition), 0);
            E_uniform_regional_alpha(i,transition,:) = sum(u.^2)*T/1001;
            E_uniform_global_alpha(i,transition) = sum(sum(u.^2))*T/1001; % integrate over inputs
        end
        
    end
    
    save(fullfile([basedir,'results/simulations/simulated_sub_Alphaeffectconc_UNIFORM_E_rho_',num2str(rho),note,'.mat']),'E_uniform_global_alpha','E_uniform_regional_alpha');

end


%%

mean_dmt = nanmean(global_CE_dmt);
mean_pcb = nanmean(global_CE_pcb);

i=1;
for rho=[1 10 20 30 40 50 60 70]
    
    load(fullfile([basedir,'results/simulations/simulated_sub_Alphaeffectconc_UNIFORM_E_rho_',num2str(rho),note,'.mat']))
    
    mean_plasma_sim = nanmean(E_uniform_global_alpha);
    
    r_dmt(i) = sqrt(sum((mean_dmt-mean_plasma_sim).^2));
    r_pcb(i) = sqrt(sum((mean_pcb-mean_plasma_sim).^2));

    i=i+1;
end

r_dmt
r_pcb

figure;
hold on
    plot(r_dmt)
    plot(r_pcb)
xlim([0.75 8.25])
xticks(1:8)
xticklabels([{'1'},{'10'},{'20'},{'30'},{'40'},{'50'},{'60'},{'70'}]);
xlabel('scaling parameter')
ylabel('Euclidean Distance')
title('Euclidean distance: AlphaConc UNIFORM')
legend('DMT vs Simulation','PCB vs Simulation')

