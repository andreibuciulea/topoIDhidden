%% Load Saved Results
%Graphical example of GSR algorithms 
open('candidate.fig');
h = gcf;
axes = get(h,'Children');
dataObjs = get(axes,'Children');

img1 = dataObjs{8}.CData;  %13,15,17,19
img2 = dataObjs{19}.CData;
img3 = dataObjs{15}.CData;
img4 = dataObjs{4}.CData;
img5 = dataObjs{6}.CData;
img6 = dataObjs{17}.CData;
img7 = dataObjs{13}.CData;
img8 = dataObjs{2}.CData;

figure()
subplot(241)
imagesc(img1/max(max(img1)))
xlabel('$ \bf{(a)}$ $\bf{S_{O}}$','Interpreter','latex','FontSize',25)
ax = gca;
ax.FontSize = 17;
subplot(242)
imagesc(img2/max(max(img2)))
xlabel(' \bf{ (b) GSt} $\bf{\hat{S}_{O}}$','Interpreter','latex','FontSize',25)
ax = gca;
ax.FontSize = 17;
subplot(243)
imagesc(img3/max(max(img3)))
xlabel(' \bf{ (c) GSm-St-GL} $\bf{\hat{S}_{O}}$','Interpreter','latex','FontSize',17)
ax = gca;
ax.FontSize = 17;
subplot(244)
imagesc(img4/max(max(img4)))
xlabel(' \bf{ (d) GSt-Rw-Fact} $\bf{\hat{S}_{O}}$','Interpreter','latex','FontSize',17)
%title('\bf{P}','Interpreter','latex','FontSize',14)
ax = gca;
ax.FontSize = 17;
colorbar()
subplot(245)
imagesc(img5/max(max(img5)))
xlabel(' \bf{(e)} $\bf{K}$','Interpreter','latex','FontSize',17)
ax = gca;
ax.FontSize = 17;
subplot(246)
imagesc(img6/max(max(img6)))
xlabel(' \bf{ (f) GSt} $\bf{\hat{K}}$','Interpreter','latex','FontSize',17)
ax = gca;
ax.FontSize = 17;
subplot(247)
imagesc(img7/max(max(img7)))
xlabel(' \bf{ (g) GSm-St-GL} $\bf{\hat{K}}$','Interpreter','latex','FontSize',17)
ax = gca;
ax.FontSize = 17;
subplot(248)
imagesc(img8/max(max(img8)))
xlabel(' \bf{ (h) GSt-Rw-Fact} $\bf{\hat{K}}$','Interpreter','latex','FontSize',17)
ax = gca;
ax.FontSize = 17;
colorbar()
sfh1 = subplot(248);
sfh1.Position = sfh1.Position + [0 0 0.001 0.001];
sfh2 = subplot(244);
sfh2.Position = sfh2.Position + [0 0 0.001 0.001];

%% Run Experiment

clear all;
close all;

addpath('../../stationarity/utils')
addpath('../../stationarity/optim')
load('exp4_S_L_data.mat');

N = 20;

g_type = 'ER';
rm_node = {'min'};
verbose = true;
verb = true;

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
Ctype = 'Cpoly'; 
C_prms.N = N;
C_prms.O = 19;
C_prms.Ctype = Ctype;
C_prms.Taps = Taps;
C_prms.noise = false;

max_iter = 10;
nG = 2;
H = 1;
models = {'GSt-Rw-Fact'};
Scores = zeros(nG,numel(Samples),numel(rm_node),H,numel(models));
opt_val = zeros(nG,1);
figure(1)
tic
for g = 1:nG
    disp(['Graph: ' num2str(g)]);
    %generate graph
    %[S, L] = generate_graph(g_type,prms);
    S = All_S(:,:,g);
    L = All_L(:,:,g);
    for s = 1:numel(Samples)
        C_prms.Samples = Samples(s);
        disp(['  Samples: ' num2str(Samples(s))]);
        C = generate_C(S,C_prms);
        all_prms = struct('M',100,'sigma',0,'norm_noise',true,'sampled',false);
        [~,~,C1] = generate_graph_signals('FA',L,all_prms,false);
        C = C/norm(C,'fro');
        for k = 1:numel(rm_node)
            links_type = rm_node{k};
            disp(['    Node selection: ' links_type]);
            for h = 1:H
                disp(['      Hidden: ' num2str(h)]);
                O = N-h;
                %select hidden nodes
                [s_n, s_h] = select_hidden_nodes(links_type, O, L, C);
                %get observed A,C,X
                Co = C(s_n,s_n);
                So = S(s_n,s_n);
                C1_o = C1(s_n,s_n);
                for m = 1:numel(models)
                    model = models{m};
                    disp(['        Model: ' model]);
                    %estimate graph
%                     reg.H = h;
                    reg = get_reg(g_type,Ctype,model,max_iter);
                    reg.S0 = So;
                    reg.H = h;
                    [S_hat,out] = estimate_graph(model,Co,reg,verbose);
                    reg1 = struct('alpha',1,'beta',0.16,'lambda',1,'gamma1',1,'gamma2',1,'epsilon',1e-6,'H',h);
                    [L1_hat,out1] = GSm_St_GL(C1_o,reg1);
                    S1_hat = diag(diag(L1_hat))-L1_hat;
                    S1_hat_bin = mbinarize(S1_hat,2);
                    %Scores(g,s,k,h+1,m+2) = out.scr;
                    %binarize S_hat
                    S_hat_bin = mbinarize(S_hat,2);
                    %compute fscore 
                    [~,~,Scores(g,s,k,h+1,m),~,~] = graph_learning_perf_eval(So,S_hat_bin);
                    [~,~,Scores1(g),~,~] = graph_learning_perf_eval(So,S1_hat_bin);
                end
            end
        end
    end
    if verb
        subplot(251)
        plot(graph(S))
        title('Graph')
        subplot(256)
        imagesc(C)
        title('Covariance matrix')
        colorbar();
        subplot(252)
        imagesc(So)
        title('S')
        colorbar()
        subplot(257)
        imagesc(C(s_n,s_h)*S(s_h,s_n))
        title('K')
        colorbar()
        subplot(255)
        imagesc(S_hat)
        title('S GSt-Rw-Fact')
        colorbar()
        subplot(2,5,10)
        imagesc(out.F_hat*out.G_hat)
        title('K GSt-Rw-Fact')
        colorbar()
        subplot(254)
        imagesc(S1_hat)
        title('S GSm-St-GL')
        colorbar()
        subplot(2,5,9)
        imagesc(out1.K)
        title('K GSm-St-GL')
        colorbar()

        figure(g+1)
        copyobj(allchild(1),g+1);
    end
end
squeeze(Scores)
toc