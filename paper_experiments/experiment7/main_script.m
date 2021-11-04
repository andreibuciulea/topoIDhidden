%% Print results
load('exp7.mat');
mrkt = {'s-','^-','^-','s-','s-'};

figure('Position',[100,100,750,550])
aux = squeeze(sum(fsc==1)/prms.nG);
prms.models ={'GSm-LR','GSm-GL','GSm-St-GL','GSm-St-LR'};
for k = [1,2,4,3]
    plot(prms.H,aux(:,k),mrkt{k},'MarkerSize',12,'LineWidth',3) 
    hold on
end    

%title('Mean Fscore')
%xlabel('Hidden variables')
%ylim([min(min(squeeze(median(fsc)))) 1])
legend(prms.models([1,2,4,3]),'FontSize',11,'FontWeight','bold')
xlabel('(c) Number of hidden variables','FontSize',17,'FontWeight','bold')
ylabel('Recovery','FontSize',17,'FontWeight','bold')
set(gca,'xtick',1:5)
ax = gca;
ax.FontSize = 17;
grid on


%%
clear
% addpath('utils/func');
% addpath('utils/sgwt_toolbox')
% addpath('utils/toolbox');
addpath(genpath('../../utils'));
addpath(genpath('../../optim'));
rng(10);
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
prms.sampled = false;
prms.L_bin = true;  % ------> Enlaces binarios

prms.nG = 128;
prms.H = 1:5; %vector de ocultas
prms.g_type = 'RBF';
prms.sig_type = 'FA';
prms.rm_node = {'rand'};
prms.models = {'GSm-LR','GSm-GL','GSm-St-LR','GSm-St-GL'};
%prms.models = {'GSm-St-LR','GSm-St-GL'};
prms.iters = 10;
prms.verb = false;

tic
[fsc,out] = get_results(prms);
toc

%%
mrkt = {'s-','^-','^-','s-','s-'};
th = 1;
figure('Position',[100,100,750,550])
aux = squeeze(sum(fsc==th)/prms.nG);
for k = 1:numel(prms.models) 
    plot(prms.H,aux(:,k),mrkt{k},'MarkerSize',12,'LineWidth',3) 
    hold on
end    

%title('Mean Fscore')
%xlabel('Hidden variables')
%ylim([min(min(squeeze(median(fsc)))) 1])
legend(prms.models,'FontSize',11,'FontWeight','bold')
xlabel('(c) Number of hidden variables','FontSize',17,'FontWeight','bold')
ylabel('Recovery','FontSize',17,'FontWeight','bold')
set(gca,'xtick',1:5)
ax = gca;
ax.FontSize = 17;
grid on


