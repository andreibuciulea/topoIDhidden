%% Plot saved results
load('data_exp5_C_Poly.mat');
th = 1;
aux = squeeze(Scores);
scores_C_poly = squeeze(sum(squeeze(Scores) ==th,1))/nG;
scores_C_poly(1,:) = [];
load('data_exp5_C_mrf.mat');
scores_C_MRF = squeeze(sum(squeeze(Scores) ==th,1))/nG;
x_data = 1:H;
mrkt = {'s--','s-','*--','*-','^--','^-','s-','s-'};
models = {'GSt-nH','GSt-Rw-nH','GSt','GSt-Rw','GSt-Rw-Fact','LVGL'};
i = 1;
figure('Position',[100,100,750,550])
for m = [3,5,6]
    plot(x_data,scores_C_poly(:,m),mrkt{i},'MarkerSize',12,'LineWidth',3);
    lgnd{i} = [replace(models{m},' ','-') ', C_{poly}'];
    i = i+1;
    hold on
    plot(x_data,scores_C_MRF(:,m),mrkt{i},'MarkerSize',12,'LineWidth',3);
    lgnd{i} = [replace(models{m},' ','-') ', C_{MRF}'];
    i = i+1;
    hold on
end
legend(lgnd,'FontSize',11,'FontWeight','bold')
xlabel('(a) Number of hidden variables','FontSize',17,'FontWeight','bold')
ylabel('Ratio of recovered graphs','FontSize',17,'FontWeight','bold')
set(gca,'xtick',1:5)
ax = gca;
ax.FontSize = 17;
grid on

%% Run Script
clear all;
close all;
%addpath('utils')
%addpath('optim')
addpath('../../../stationarity/utils');
addpath('../../../stationarity/optim');
%rng(23)

N = 20;
Obs = 19;

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
Samples = 1e6;
Ctypes = {'Cmrf','Cpoly'};
noise = false;
max_iter = 5;
nG = 32;
H = 5;

models = {'GSt','GSt-Rw-Fact','LVGL'};
Scores = zeros(nG,numel(Ctypes),H,numel(models));

opt_val = zeros(nG,1);
tic
parfor g = 1:nG
    disp(['Graph: ' num2str(g)]);
    %generate graph
    [S, L] = generate_graph(g_type,prms);
    Scores_g = zeros(numel(Ctypes),H,numel(models));
    C_prms = struct('N',N,'Taps',Taps,'noise',noise,'Samples',Samples); 
    for ct = 1:numel(Ctypes)
        C_prms.Ctype = Ctypes{ct};
        C = generate_C(S,C_prms);
        for h = 1:H
            %disp(['      Hidden: ' num2str(h)]);
            O = N-h;
            %select hidden nodes
            [s_n, s_h] = select_hidden_nodes(rm_node, O, L, C);
            %get observed A,C,X
            Co = C(s_n,s_n);
            So = S(s_n,s_n);
            for m = 1:numel(models)
                model = models{m};
                %disp(['        Model: ' model]);
                %estimate graph
                reg = get_reg(g_type,Ctypes{ct},model,max_iter);
                reg.H = h;
                [S_hat,out] = estimate_graph(model,Co,reg,verbose);
                %Scores(g,s,k,h+1,m+1) = out.scr;
                %binarize S_hat
                S_hat_bin = mbinarize(S_hat,2);
                %compute fscore 
                [~,~,Scores_g(ct,h,m),~,~] = graph_learning_perf_eval(So,S_hat_bin);   
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
x_data = 1:H;
mrkt = {'s--','s-','*--','*-','^--','^-','s-','s-'};
%models = {'GSt-nH','GSt-Rw-nH','GSt','GSt-Rw','GSt-Rw-Fact','LVGL'};
i = 1;
figure('Position',[100,100,750,550])
for m = 1:numel(models)
    plot(x_data,scores_C_poly(:,m),mrkt{i},'MarkerSize',12,'LineWidth',3);
    lgnd{i} = [replace(models{m},' ','-') ', C_{poly}'];
    i = i+1;
    hold on
    plot(x_data,scores_C_MRF(:,m),mrkt{i},'MarkerSize',12,'LineWidth',3);
    lgnd{i} = [replace(models{m},' ','-') ', C_{MRF}'];
    i = i+1;
    hold on
end
legend(lgnd,'FontSize',11,'FontWeight','bold')
xlabel('(a) Number of hidden variables','FontSize',17,'FontWeight','bold')
ylabel('Ratio of recovered graphs','FontSize',17,'FontWeight','bold')
set(gca,'xtick',1:5)
ax = gca;
ax.FontSize = 17;
grid on
