
N = 20;
links_type = 'rand';
g_type = 'ER';
verb = false;

Taps = 5;
Samples = 1e6;
Ctype = 'poly'; 
C_prms.N = N;
C_prms.Samples = Samples;
C_prms.Ctype = Ctype;
C_prms.Taps = Taps;
C_prms.noise = false;

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

nG = 12;

epsils1 = 1e-6;  %1e-6
epsils2 = [1];  %1e-6
dels = [1];  %1e-6
taus = [1];  %1e-6
n_e1 = numel(epsils1);
n_e2 = numel(epsils2);
n_d = numel(dels);
n_t = numel(taus);

O = 19;
C_prms.O = O; %-----------------------> cuidado con esto
max_iter = 10;
model = 'GSR v1';

Scores = zeros(nG,n_e1,n_e2,n_d,n_t);
tic
for g = 1:nG
    disp(['Graph: ' num2str(g)]);
    [S, L] = generate_graph(g_type,prms);
    [C] = generate_C(S,C_prms);
    [s_n, s_h] = select_hidden_nodes(links_type, O, L, C);
    Co = C(s_n,s_n);
    So = S(s_n,s_n);
    Scores_g = zeros(n_e1,n_e2,n_d,n_t);
    for e1 = 1:n_e1
        for e2 = 1:n_e2
            for d = 1:n_d
                for t = 1:n_t
                    reg = struct('eps1',epsils1(e1),'eps2',epsils2(e2),'del',...
                                  dels(d),'tau',taus(t),'max_iter',max_iter);
                    [S_hat, ~] = estimate_graph(model,Co,reg,verb);
                    S_hat_bin = mbinarize(S_hat,2);
                    [~,~,Scores_g(e1,e2,d,t),~,~] = graph_learning_perf_eval(So,S_hat_bin);
                end
            end
        end
    end
    Scores(g,:,:,:,:) = Scores_g;
    
end
disp(['Elapsed time: ', num2str(toc/60) ' min'])


%%
figure()
subplot(311)
imagesc(squeeze(mean(Scores)))
set(gca,'XTick',[1:numel(epsils2)], 'XTickLabel', epsils2);
set(gca,'YTick',[1:numel(epsils1)], 'YTickLabel', epsils1);
set(gca,'XTickLabelRotation',45);
xlabel('epsils2')
ylabel('epsils1')
colorbar()
subplot(312)
imagesc(squeeze(median(Scores)))
colorbar()
subplot(313)
imagesc(squeeze(mean(Scores==1)))
colorbar()

%%
epsils1 = logspace(0,-6,7);
data1 = {'gsr_v1_fro.mat','gsr_v1_abs.mat','gsr_v1_nada.mat'};
data2 = {'gsr_v1_fro_1e-2.mat','gsr_v1_abs_1e-2.mat','gsr_v1_nada_1e-2.mat'};
figure()
for k = 1:3
    load(data1{k})
    semilogx(epsils1,median(Scores))
    hold on
end
legend('fro','abs1','nada1');

figure()
for k = 1:3
    load(data2{k})
    semilogx(epsils1,median(Scores))
    hold on
end
legend('fro2','abs2','nada2');


