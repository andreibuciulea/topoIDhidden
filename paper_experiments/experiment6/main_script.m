%% Plot saved results
load('poly_exp6.mat');
th = 1;
aux = squeeze(Scores);
scores_C_poly = squeeze(sum(squeeze(Scores(:,:,:,2,:)) == th,1))/nG;
%scores_C_poly(1,:) = [];
load('cmrf_exp6.mat');
scores_C_MRF = squeeze(sum(squeeze(Scores(:,:,:,2,:)) == th,1))/nG;
x_data = Samples;
mrkt = {'s--','s-','*--','*-','^--','^-','s-','s-'};
i = 1;
models = {'','','GSt','','GSt-Rw-Fact','LV-GL'};
figure('Position',[100,100,750,550])
for m = [3,5,6]
    semilogx(x_data,scores_C_poly(:,m),mrkt{i},'MarkerSize',12,'LineWidth',3);
    lgnd{i} = [replace(models{m},' ','-') ', C_{poly}'];
    i = i+1;
    hold on
    semilogx(x_data,scores_C_MRF(:,m),mrkt{i},'MarkerSize',12,'LineWidth',3);
    lgnd{i} = [replace(models{m},' ','-') ', C_{MRF}'];
    i = i+1;
    hold on
end
lg = legend(lgnd,'FontSize',11,'FontWeight','bold');
set(lg,'color','none');
xlabel('(b) Number of samples','FontSize',17,'FontWeight','bold')
ylabel('Recovery','FontSize',17,'FontWeight','bold')
set(gca,'xtick',Samples)
grid on
ax = gca;
ax.FontSize = 17;
xticks([1e2,1e3,1e4,1e5,1e6]);
xticklabels({'10^2','10^3','10^4','10^5','10^6'});

%%
clear all;
close all;
addpath('../../../stationarity/utils');
addpath('../../../stationarity/optim');
%rng(23)

N = 20;

g_type = 'ER';
rm_node = 'min';
verbose = false;

prms.N = N;
prms.M = 1;
prms.p = 0.2;%ER
prms.m = 1;%BA
prms.m0 = 1;%BA
prms.K = 4; %SW
prms.beta = 0.3;%SW
prms.connected = true; %RBF
prms.norm_noise = true;
prms.norm_L = false;
prms.sigma = 0;

Taps = 5;
Samples = round(logspace(2,6,10));
Ctypes = {'Cmrf','Cpoly'}; 
noise = true;

max_iter = 5;
nG = 32;
H = 1;
O = N-H;
models = {'GSt','GSt-Rw-Fact','LVGL'};
Scores = zeros(nG,numel(Ctypes),numel(Samples),numel(models));
opt_val = zeros(nG,1);
tic
parfor g = 1:nG
    disp(['Graph: ' num2str(g)]);
    %generate graph
    [S, L] = generate_graph(g_type,prms);
    Scores_g = zeros(numel(Ctypes),numel(Samples),numel(models));
    C_prms = struct('N',N,'Taps',Taps,'noise',noise);
    for s = 1:numel(Samples)
        smpl = Samples(s);
        C_prms.Samples = smpl;
        %disp(['  Samples: ' num2str(smpl)]);
        for ct = 1:numel(Ctypes)
            C_prms.Ctype = Ctypes{ct};
            C = generate_C(S,C_prms);
            %select hidden nodes
            [s_n, s_h] = select_hidden_nodes(rm_node, O, L, C);
            %get observed A,C,X
            Co = C(s_n,s_n);
            So = S(s_n,s_n);
            for m = 1:numel(models)
                model = models{m};
                %disp(['        Model: ' model]);
                reg = get_reg(g_type,Ctypes{ct},model,max_iter);
                %reg.S0 = So;
                reg.H = H;
                [S_hat,out] = estimate_graph(model,Co,reg,verbose);
                %Scores(g,s,k,h+1,m+1) = out.scr;
                %binarize S_hat
                S_hat_bin = mbinarize(S_hat,2);
                %compute fscore 
                [~,~,Scores_g(ct,s,m),~,~] = graph_learning_perf_eval(So,S_hat_bin);
            end
        end
    end
    Scores(g,:,:,:) = Scores_g;
end
disp([num2str(toc/60) 'min']);
%%
th = 1;
scores_C_MRF = squeeze(sum(squeeze(Scores(:,1,:,:))==th,1))/nG;
scores_C_poly = squeeze(sum(squeeze(Scores(:,2,:,:))==th,1))/nG;
x_data = Samples;
mrkt = {'s--','s-','*--','*-','^--','^-','s-','s-'};
i = 1;
%models = {'','','GSt','','GSt-Rw-Fact','LV-GL'};
figure('Position',[100,100,750,550])
for m = 1:numel(models)
    semilogx(x_data,scores_C_poly(:,m),mrkt{i},'MarkerSize',12,'LineWidth',3);
    lgnd{i} = [replace(models{m},' ','-') ', C_{poly}'];
    i = i+1;
    hold on
    semilogx(x_data,scores_C_MRF(:,m),mrkt{i},'MarkerSize',12,'LineWidth',3);
    lgnd{i} = [replace(models{m},' ','-') ', C_{MRF}'];
    i = i+1;
    hold on
end
lg = legend(lgnd,'FontSize',11,'FontWeight','bold');
set(lg,'color','none');
xlabel('(b) Number of samples','FontSize',17,'FontWeight','bold')
ylabel('Recovery','FontSize',17,'FontWeight','bold')
set(gca,'xtick',Samples)
grid on
ax = gca;
ax.FontSize = 17;
xticks([1e2,1e3,1e4,1e5,1e6]);
xticklabels({'10^2','10^3','10^4','10^5','10^6'});

