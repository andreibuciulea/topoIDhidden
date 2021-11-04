%% Load Saved Results
load('exp1.mat');
fmts = {'s-','x-','o-','*-','o','o-'};
figure('Position',[100,100,750,550])
mdls = {'GL-SigRep','GSm','GSm-LR','GSm-GL'};
for k = [1,2,4,6]
    plot(prms.H,squeeze(median(fsc(:,:,:,k))),fmts{k},'MarkerSize',12,'LineWidth',3)
    hold on
end
legend(mdls,'FontSize',11,'FontWeight','bold')
%title('Median Fscore')
xlabel('(a) Number of hidden variables','FontSize',17,'FontWeight','bold')
ylabel('Fscore','FontSize',17,'FontWeight','bold')
set(gca,'xtick',1:5)
ylim([0.89 0.92])
ax = gca;
ax.FontSize = 17;
grid on

%% Run Experiment
clear all
rng(10)
addpath(genpath('../../utils'));
addpath(genpath('../../optim'));

prms.N = 20;
prms.norm_L = true;
prms.sigma = 0;
prms.p = 0.3;%ER
prms.m = 1;%BA
prms.connected = true; %RBF
prms.K = 4; %SW
prms.beta = 0.3;%SW
prms.M = 100;
prms.norm_noise = true;
prms.sampled = true;
prms.L_bin = true; 

prms.nG = 1024;
prms.H = 1:5;
prms.g_type = 'RBF';
prms.sig_type = 'FA';
prms.rm_node = {'rand'};
prms.models = {'GL-SigRep','GSm','GSm-LR','GSm-GL'};
prms.iters = 1;
prms.verb = false;

tic
[fsc,out] = get_results(prms);
toc

%% Print Results
fmts = {'s-','x-','*-','o-'};
figure('Position',[100,100,750,550])
mdls = prms.models;
for k= 1:numel(mdls)
    plot(prms.H,squeeze(median(fsc(:,:,:,k))),fmts{k},'MarkerSize',12,'LineWidth',3)
    hold on
end
legend(mdls,'FontSize',11,'FontWeight','bold')
xlabel('(a) Number of hidden variables','FontSize',17,'FontWeight','bold')
ylabel('Fscore','FontSize',17,'FontWeight','bold')
set(gca,'xtick',1:5)
ax = gca;
ax.FontSize = 17;
grid on

