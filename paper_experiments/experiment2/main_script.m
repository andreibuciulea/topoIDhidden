%% ER graphs with p = [0.1,0.3,0.5] varying the noise level.
%% Load Saved Results
load('exp2.mat');
fmts = {'s-','v-','s--','v--','s:','v:','*-','*--','+-','+--','o:','s:','x:','*:','+:','<:'};
trh = 0.9;
mdl = {'GSm-LR','GSm-GL'};
i = 1;
%all_fsc = out;
figure('Position',[100,100,750,550])
for gp = [1,2,3] %probabilities
    %figure()
    for k = [1,2] %models
        plot(sigmas,squeeze(median(all_fsc(gp,:,:,:,k),4)),fmts{i},'MarkerSize',12,'LineWidth',3)
        hold on
        lgnd{i} = [mdl{k} ', ER p=' num2str(probs(gp))];
        i = i+1;
    end
end
legend(lgnd,'FontSize',11,'FontWeight','bold')
%title(['Median Fscore ' graphs_type{gt} ' graph'])
xlabel('(b) Normalized noise power','FontSize',17,'FontWeight','bold')
ylabel('Fscore','FontSize',17,'FontWeight','bold')
ax = gca;
ax.FontSize = 17;
grid on
ylim([0.60 0.97])

%% Run Experiment
clear all

addpath(genpath('../../utils'));
addpath(genpath('../../optim'));

prms.N = 20;
prms.norm_L = true;
prms.sigma = 0;
prms.p = 0.2;%ER
prms.m = 1;%BA
prms.connected = true; %RBF
prms.K = 3; %SW
prms.beta = 0.3;%SW
prms.M = 100;
prms.norm_noise = true;
prms.sampled = true;
prms.L_bin = true;

prms.nG = 32;%128;
prms.H = 1:1; %
prms.sig_type = 'FA';
prms.rm_node = {'rand'};
prms.models = {'GSm-LR','GSm-GL'};
prms.iters = 1;
prms.verb = false;

graphs_type = {'ER'};
sigmas = [0:0.02:0.2];
probs = [0.1,0.3,0.5];
all_fsc = zeros(numel(probs),numel(graphs_type),numel(sigmas),prms.nG,numel(prms.models));
all_out = zeros(numel(probs),numel(graphs_type),numel(sigmas),prms.nG,numel(prms.models));

tic
for sp = 1:numel(probs)
    prms.p = probs(sp);
    for gt = 1:numel(graphs_type)
        prms.g_type = graphs_type{gt};
        for sg = 1:numel(sigmas)
            prms.sigma = sigmas(sg);
            [fsc,out] = get_results(prms);
            all_fsc(sp,gt,sg,:,:) = squeeze(fsc);
            %all_out(sp,gt,sg,:,:) = squeeze(out);
        end
    end
end
toc

%% Print Results
fmts = {'s-','v-','s--','v--','s:','v:','*-','*--','+-','+--','o:','s:','x:','*:','+:','<:'};
mdl = prms.models;
i = 1;
figure()

for pr = 1:numel(probs)
    for k = 1:numel(mdl)
        plot(sigmas,squeeze(median(all_fsc(pr,:,:,:,k),4)),fmts{i},'MarkerSize',10,'LineWidth',2)
        hold on
        lgnd{i} = [mdl{k} ', ER p=' num2str(probs(pr))];
        i = i+1;
    end
    legend(lgnd,'FontSize',11,'FontWeight','bold')
    xlabel('(b) Normalized noise power','FontSize',17,'FontWeight','bold')
    ylabel('Fscore','FontSize',17,'FontWeight','bold')
    ax = gca;
    ax.FontSize = 17;
end

