load('st_data/st_data.mat');
load('st_data\ave_real_data_prediction.mat')
%%
nG = 15;
A0 = contact_matrices{nG};
C0 = mi_matrices{nG};
sg = 7:47;
A = A0(sg,sg);
C = C0(sg,sg);
%C2 = C0(7:48,7:48);
%C3 = C0(9:46,9:46);

N = size(A,1);
M = 250;

%Dong --> 
%models = {'adj-eig','adj-C','adj-C-H','estim-LC','hLAo','LnhC-st-log-nuc','LV GL'};
models = {'LnhC-st-log-nuc'};
results = zeros(M,numel(models));
verbose = false;
max_iter = 5;
H = 12; %poner H lo que sea 4,5,6
figure()
for k = 1:numel(models)
    model = models{k};
    switch model
        case 'estim-LC'
            reg = struct('alpha',0.1,'beta',0.1);
            L_hat = estimate_L_C(C,reg);
            A_hat = diag(diag(L_hat))-L_hat;
            subplot(221)
            imagesc(A_hat)
            title(model)
        case 'hLAo'
            %reg = struct('alpha',1,'beta',0.01,'lambda',0.1,'gamma1',0.1);
            reg = struct('alpha',1,'beta',1e-3,'lambda',1e-3,'gamma1',1,'epsilon',10);
            [L_hat, ~, ~] = hsmooth_LAo(C,H,reg,verbose);
            A_hat = diag(diag(L_hat))-L_hat;
%             subplot(223)
%             imagesc(A_hat)
%             title(model)
        case 'GSm-St-LR'
            reg = struct('alpha',1,'beta',0.081,'lambda',0.658,'gamma1',10,'epsilon',0.8,'H',12);
            [L_hat,~] = estimate_L_C_st_log_nuc(C,reg);
            A_hat = diag(diag(L_hat))-L_hat;
        case 'GSm-St-GR'
            reg = struct('alpha',1,'beta',1,'lambda',1,'gamma1',1,'gamma2',1,'epsilon',0.7,'H',12);
            [L_hat,~] = estimate_L_C_st_log_nuc(C,reg);
            A_hat = diag(diag(L_hat))-L_hat;
        case 'adj-eig'
            reg = struct('epsil_min', 0.95,'epsil_max',0.95,...
                'alpha',1,'beta',1,'max_iter',max_iter);
            A_hat = our_recovery_adjacency(C, reg);
        case 'adj-C'
            reg = struct('epsil_min', 0.4,'epsil_max',0.4,...
                'alpha',1,'beta',1,'max_iter',max_iter);
            A_hat = our_recovery_adjacency_v2(C, reg);
        case 'adj-C-H'
            reg = struct('epsil_min', 0.9,'epsil_max',0.9,...
                'alpha',1.269e-5,'beta',1,'max_iter',max_iter);
            [A_hat,~] = our_recovery_adjacency_v3(C, reg);
            subplot(222)
            imagesc(A_hat)
            title(model)
        case 'adj-GSt'
            reg = struct('alpha',6.8,'beta',3.5,'epsilon',0.8,'epsil_min', 0.8,'epsil_max',0.8);
            [A_hat,~]= our_recovery_GSt(C, reg)
        case 'GSR-Rw-Fact'
            reg = struct('alpha',8.11e3,'beta',1.52e4,'tau',1,'del',9e-3,'R',11,'max_iter',5,'epsilon',0.1);
            [A_hat,~] = GSR_H_Rw_Fact_real_data(C,reg,verbose);
        case 'LV GL'
            reg = struct('alpha',1e-4,'beta',1e-4);
            [A_hat, ~] = LVGLASSO(C,reg,verbose);
            subplot(224)
            imagesc(A_hat)
            title(model)
        otherwise
            disp('Unknown optimization model');
            L_hat = NaN;
    end
    results(:,k) = obtain_est_results(A_hat,A);
end

%
figure()
for k = 1:numel(models)
    plot(results(:,k),'LineWidth',2);
    hold on
end
legend(models)
figure()
imagesc(A_hat)
%%
plot(mi_ave)
hold on
plot(di_ave)
hold on
plot(nd_ave)
hold on
plot(nd_di_ave)
legend([models,'MI','DI','ND','DI+ND'])


%% Check Params
load('../real_data/st_data.mat');
nG = 15;
A0 = contact_matrices{nG};
C0 = mi_matrices{nG};
sg = 7:47;
A = A0(sg,sg);
C = C0(sg,sg);

alphas = 1;
betas = logspace(0,-6,2);
lambdas = logspace(0,-6,2);
gammas1 = logspace(0,-6,2);
epsilon = 0.8;

n_a = numel(alphas);
n_b = numel(betas);
n_la = numel(lambdas);
n_ga = numel(gammas1);

models = {'GSm-GL'};
n_mo = numel(models);
results = zeros(n_b,n_a,n_la,n_ga,n_mo);


max_iter = 5;
H = 1;
verbose = false;

tic
for b = 1:n_b
    results_b = zeros(n_a,n_la,n_ga,n_mo);
    A_hat = zeros(size(A));
    for a = 1:n_a
        for la = 1:n_la
            for ga = 1:n_ga
                for mo = 1:n_mo
                    model = models{mo};
                    reg = struct('alpha',alphas(a),'beta',betas(b),'lambda',lambdas(la),'gamma1',gammas1(ga),...
                        'epsil_min',epsilon,'epsil_max',epsilon,'epsilon',epsilon,'max_iter',max_iter,...
                        'tau',1,'del',1e-3,'H',H);
                    switch model
                        case 'estim-LC'
                            L_hat = estimate_L_C(C,reg);
                            A_hat = diag(diag(L_hat))-L_hat;
                        case 'hLAo'
                            [L_hat, ~, ~] = hsmooth_LAo(C,H,reg,verbose);
                            A_hat = diag(diag(L_hat))-L_hat;
                       case 'GSm-GL'
                            [L_hat, ~, ~] = GSm_GL(C,H,reg,false);
                        case 'GSm-St-LR'
                            [L_hat,~] = GSm_St_LR(C,reg);
                            A_hat = diag(diag(L_hat))-L_hat;
                        case 'GSm-St-GL'
                            [L_hat,~] = GSm_St_GL(C,reg);
                            A_hat = diag(diag(L_hat))-L_hat;
                        case 'adj-eig'
                            A_hat = our_recovery_adjacency(C, reg);
                        case 'adj-C'
                            A_hat = our_recovery_adjacency_v2(C, reg);
                        case 'adj-GSt'
                            [A_hat,~]= our_recovery_GSt(C, reg);
                        case 'adj-C-H'
                            A_hat = our_recovery_adjacency_v3(C, reg);
                        case 'LV GL'
                            [A_hat, ~] = LVGLASSO(C,reg,verbose);
                        case 'GSR-Rw-Fact'
                            %reg = struct('alpha',100,'beta',1,'tau',1,'del',1e-3,'R',11,'max_iter',5);
                            [A_hat,~] = GSR_H_Rw_Fact_real_data(C,reg,verbose);
                        otherwise
                            disp('Unknown optimization model');
                            L_hat = NaN;
                    end
                    res = obtain_est_results(A_hat,A);
                    results_b(a,la,ga,ep) = res(end); 
                end
            end
        end
    end

    results(b,:,:,:,:) = results_b;
end

disp(['Elapsed time: ' num2str(toc/60) 'min'])



