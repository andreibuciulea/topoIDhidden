%%Load results
load('exp3.mat')
figure('Position',[100,100,750,550])
plot(smth,fsc,'*')
xlabel('(c) Local Variation (LV)','FontSize',17,'FontWeight','bold')
ylabel('Fscore','FontSize',17,'FontWeight','bold')
xlim([0 2.2])
ax = gca;
ax.FontSize = 17;
grid on


%% Signals with different total variation values
clear all

addpath(genpath('../../utils'));
addpath(genpath('../../optim'));

R = 128;%128;
N = 30;  %30
M = 1e4; %1e4
prms.nG = 1;
h = 1;
prms.M = M;
prms.N = N;
prms.norm_L = false;
prms.L_bin = true;
prms.connected = true;
g_type= 'RBF';
sig_type = 'FA1';
prms.sigma = 0;
prms.sampled = true;
prms.g_type = g_type;
prms.sig_type = sig_type;
rm_node = 'rand';
model = 'GSm-LR';
iters = 10;
verb = false;

[A,L] = generate_graph(g_type,prms);
[V, Lambda] = eig(L);

rng(1)

smth = zeros(R,N);
smth_o = zeros(R,N);
fsc = zeros(R,N);
tic
parfor r = 1:R
    %[X,~,C,~] = generate_graph_signals(sig_type, L, prms, false);
    h_mu = zeros(N,1);
    fsc_r = zeros(N,1);
    smth_r = zeros(N,1);
    smth_o_r = zeros(N,1);
    for k = 3:N-2
        H = mvnrnd(h_mu, eye(N), M)';
        V2 = V(:,k-2:k+2); H2 = H(k-2:k+2,:);%./vecnorm(H(k,:));
        X = V2*H2;  C = X*X'/M;
        %%%%%
        
        smth_r(k) = vec(C)'*vec(L)/N; 
        %select hidden node
        [s_n, s_h] = select_hidden_nodes(rm_node, N-h, L, X);
        Xo = X(s_n,:); Ao = A(s_n,s_n); Co = C(s_n,s_n); Lo = L(s_n,s_n);
        smth_o_r(k) = vec(Co)'*vec(Lo)/N; 
        % get regularizers
        reg = get_reg(g_type,sig_type,model);
        [L_hat, ~] = estimate_graph(model, Xo, Co, h, reg, iters, verb);
        A_aux = diag(diag(L_hat))-L_hat;
        A_hat = mbinarize(A_aux,2);
        Ao = mbinarize(Ao,2);
        [~,~,fsc_r(k),~,~] = graph_learning_perf_eval(Ao,A_hat);
    end
    fsc(r,:) = fsc_r;
    smth(r,:) = smth_r;
    smth_o(r,:) = smth_o_r;
end
toc
%%
figure('Position',[100,100,750,550])
plot(smth,fsc,'*')
xlabel('(c) Local Variation (LV)','FontSize',17,'FontWeight','bold')
ylabel('Fscore','FontSize',17,'FontWeight','bold')
xlim([0 2.2])
ax = gca;
ax.FontSize = 17;
grid on
