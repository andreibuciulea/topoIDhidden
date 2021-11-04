rng(1)
N = 20;
links_type = 'min';
g_type = 'ER';
verb = false;

Taps = 5;
Samples = 1e6;
Ctype = 'Cpoly'; 
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

nG = 32;
%ajustar el GSR H v1 para fro y resta y ver si hay diferencia
alphas = logspace(2,-2,10);  %1e-6
betas = logspace(3,1,6); %1e-6
epsils1 = 1e-6;  %1e-6
epsils2 = 5e-2;  %1e-6
tau = 1;
del = 1e-3;

n_a = numel(alphas);
n_b = numel(betas);
n_e1 = numel(epsils1);
n_e2 = numel(epsils2);

O = 19;
C_prms.O = O; 
max_iter = 5;
model = 'GSt-Rw-Fact';

Scores = zeros(nG,n_a,n_b,n_e1,n_e2);
tic
parfor g = 1:nG
    disp(['Graph: ' num2str(g)]);
    [S, L] = generate_graph(g_type,prms);
    [C] = generate_C(S,C_prms);
    
    %disp('Normalized C')
    %C = C/norm(C,'fro');
    [s_n, s_h] = select_hidden_nodes(links_type, O, L, C);
    Co = C(s_n,s_n);
    So = S(s_n,s_n);
    Scores_g = zeros(n_a,n_b,n_e1,n_e2);
    for a = 1:n_a
        for b = 1:n_b
            for e1 = 1:n_e1
                for e2 = 1:n_e2
                    reg = struct('eps1',epsils1(e1),'eps2',epsils2(e2),'alpha',alphas(a),...
                                 'beta',betas(b),'max_iter',max_iter,'g_type', g_type,...
                                 'Ctype', Ctype,'tau',tau,'del',del,'H',N-O,'S0',So);
                    [S_hat, ~] = estimate_graph(model,Co,reg,verb);
                    S_hat_bin = mbinarize(S_hat,2);
                    [~,~,Scores_g(a,b,e1,e2),~,~] = graph_learning_perf_eval(So,S_hat_bin);
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
imagesc(squeeze(sum(Scores==1)))
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


