%% Print results

for k = 1:numel(models)
    disp(['Model: ', models{k},' Fscore: ', num2str(fsc(k)),' Precision: ', num2str(prec(k)),' Recall: ', num2str(rec(k)),' NMI: ', num2str(nmi(k))]);
end

%% Run experiment
% clear all
close all
filename = 'real_data_meteo.xlsx';
%data_alt = xlsread(filename);
table = readtable(filename);
data_alt = table2array(table(:,3));
signals = str2double(table2array(table(:,7:18)));

max_alt = 300;
N = numel(data_alt);
S88 = zeros(N,N);

for i = 1:N
    for j = 1:N
        if abs(data_alt(i)-data_alt(j)) < max_alt
            S88(i,j) = 1;
        end
    end
end
S88 = S88-diag(diag(S88));
G88 = graph(S88);

N = 20;
S = zeros(N,N);
for i = 1:N
    for j = 1:N
        if abs(data_alt(i)-data_alt(j)) < max_alt
            S(i,j) = 1;
        end
    end
end

S = S-diag(diag(S));
G20 = graph(S);
% figure()
% subplot(2,2,1)
% imagesc(S88)
% title('Adj matrix 88 stations')
% subplot(2,2,2)
% plot(G88)
% title('Graph 88 nodes')
% subplot(2,2,3)
% imagesc(S)
% title('Adj matrix 20 stations')
% subplot(2,2,4)
% plot(G20)
% title('Graph 20 nodes')
% 
% figure()
% imagesc(signals)
% title('Tempetarure of 88 stations')
% xlabel('Months')
% ylabel('Strations')
% colorbar()


%%
models = {'GL-SigRep','GSm-LR','GSm-GL','GSm-St-LR','GSm-St-GL'};
N = 20;
signals_N = signals(1:N,:);
%signals_N = signals_N - mean(signals_N);
%signals_N = signals_N/max(max(signals_N));
X = signals_N;
L_N = diag(sum(S))-S;
C = signals_N*signals_N';
%
C = C/norm(C,'fro');
epsilons = 0.9;
prec = zeros(numel(models),1);
rec = zeros(numel(models),1);
fsc = zeros(numel(models),1);
nmi = zeros(numel(models),1);

for eps = 1:numel(epsilons)
    for k = 1:numel(models)
        model = models{k};
        reg = real_data_reg(model);
        reg.H = 1;
        %reg.epsilon = epsilons(eps);
        [L_hat,cond(eps,:),R(eps,:,:)] = real_data_estim(X,C,model,1,reg);
        S0 = diag(diag(L_hat))-L_hat;
        [S_hat] = mbinarize(S0,2);
        L_hat_bin = diag(sum(S_hat)) - S_hat;
        [prec(k,eps),rec(k,eps),fsc(k,eps),nmi(k,eps),~] = graph_learning_perf_eval(S,S_hat);
        disp(['Model: ', model, ' Epsilon: ', num2str(epsilons),' Fscore: ', num2str(fsc(k,eps))]);
        all_L(eps,:,:) = L_hat; 
        all_C(eps,:,:) = C;
    end
end






