clear all;
close all;
addpath('utils')
addpath('optim')
%rng(23)

N = 20;

g_type = 'ER';
rm_node = {'min'};
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
Ctype = 'poly'; 
C_prms.N = N;
C_prms.O = 20; %-----------------------> cuidado con esto (para python y regain)
C_prms.Ctype = Ctype;
C_prms.Taps = Taps;
C_prms.noise = false;

max_iter = 5;
nG = 1;
H = 1;
%models = {'GSR', 'GSR v1', 'GSR Rw', 'GSR H','GSR H v1'};
%models = {'GSR v1','GSR Rw v1','GSR H v1','GSR H Rw v1'};
%models = {'GSR','GSR Rw','GSR H','GSR H Rw','GSR H Rw Fact','GSR H Rw Full'};
models = {'GSR H Rw Fact'};
%
Scores = zeros(nG,numel(Samples),numel(rm_node),H,numel(models));
opt_val = zeros(nG,1);

for g = 1:nG
    disp(['Graph: ' num2str(g)]);
    %generate graph
    [S, L] = generate_graph(g_type,prms);
    %generate covariance (se podría poner otro bucle en función de las realizaciones)
    for s = 1:numel(Samples)
        C_prms.Samples = Samples(s);
        disp(['  Samples: ' num2str(Samples(s))]);
        C = generate_C(S,C_prms);
        %C = C/norm(C,'fro');
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
                for m = 1:numel(models)
                    model = models{m};
                    disp(['        Model: ' model]);
                    %estimate graph
%                     reg.H = h;
                    reg = get_reg(g_type,Ctype,model,max_iter);
                    reg.S0 = So;
                    reg.R = h;
                    [S_hat,out] = estimate_graph(model,Co,reg,verbose);
                    %Scores(g,s,k,h+1,m+1) = out.scr;
                    %binarize S_hat
                    S_hat_bin = mbinarize(S_hat,2);
                    %compute fscore 
                    [~,~,Scores(g,s,k,h+1,m),~,~] = graph_learning_perf_eval(So,S_hat_bin);
                    
                end
            end
        end
    end
end
squeeze(Scores)
%%
%clear plt_data
plot_type = 2; 
s = 1;
h = 1;
leg = models;
tit = rm_node;
x_lab = 'Hidden variables';
y_lab = 'Fscore';
measure_types = {'mean','median','recovery'};
mrkt = {'s:','s:','^--','^-','*-','+-','s-','s-'};
f_clr = [0    1     0     1    0    0]; 
x_data = 1:H;

for m = 1:numel(measure_types)
    switch measure_types{m}
        case 'mean'
            scores = squeeze(mean(Scores,1));
        case 'median'
            scores = squeeze(median(Scores,1));
        case 'recovery'
            aux = squeeze(Scores);
            scores = squeeze(sum(aux ==1,1))/nG;
        otherwise
            scores = -1;
    end
    scores(1,:) = [];
    if plot_type == 1
        %plots all node selection for differen models
        plt_data(:,:,1) = squeeze(scores(s,:,:))';
        tit = models;
        leg = rm_node;
    elseif plot_type == 2
        %plots all models for different node selection
        if numel(size(scores)) == 2
            plt_data(:,:,1) = scores;
        else
             plt_data(1,:,:) = squeeze(scores(s,:,:))';
        end
    elseif plot_type == 3
        plt_data(1,:,:) = squeeze(scores(:,:,h));
        x_lab = 'Samples';
        x_data = Samples;
    end
    
    figure()
    
    for t = 1:size(plt_data,2)
        if plot_type == 3
            semilogx(x_data,plt_data(:,:,t))
        else
           hold on
           fgr = plot(x_data,plt_data(:,t),mrkt{t},'MarkerSize',8,'LineWidth',1.5);
           if f_clr(t)
               fgr.MarkerFaceColor = fgr.Color;
           end
        end
    end
    xlabel(x_lab)
    ylabel(y_lab)
    title(['Ms type: ' char(measure_types{m}) ' Nd sel: ' char(tit)])
    legend(leg)  
end


