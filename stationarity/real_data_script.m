clear all;
close all;
addpath('utils')
addpath('optim')

rm_node = {'rand'};
verbose = false;

alphas = 1.269e-5;
betas = 1;
epsilons = 0.7:0.1:1;%(0.85-->v3) %(0.9 --> v2)
%en el v2 me faltara buscar el epsilon y ya está
n_a = numel(alphas);
n_b = numel(betas);
n_e = numel(epsilons);

max_iter = 5;
models = {'simple rew'};

load('st_data/st_data.mat');
nG = 15;
A0 = contact_matrices{nG};
C0 = mi_matrices{nG};
sg = 1:53;
A = A0(sg,sg);
C1 = C0(sg,sg);
C2 = C0(7:48,7:48);
C3 = C0(9:46,9:46);
% norm(C*A-A*C,'fro');
% C = C/norm(C,'fro');
% A = A/norm(A,'fro');
% norm(C*A-A*C,'fro');
N = size(A,1);
A_hat_all_1 = zeros(4,n_a*n_b*n_e,size(C1,1),size(C1,1));
A_hat_all_2 = zeros(4,n_a*n_b*n_e,size(C2,1),size(C2,1));
A_hat_all_3 = zeros(4,n_a*n_b*n_e,size(C3,1),size(C3,1));
%%
i = 1;
tic
for a = 1:n_a
    for b = 1:n_b
        for e = 1:n_e
            reg = struct('epsil_min', epsilons(e),'epsil_max',epsilons(e),'alpha',alphas(a),...
                'beta',betas(b),'max_iter',max_iter);
%             A_hat_all_1(1,i,:,:) = our_recovery_adjacency(C1, reg);
%             A_hat_all_1(2,i,:,:) = our_recovery_adjacency_v2(C1, reg);
%             A_hat_all_1(3,i,:,:) = our_recovery_adjacency_v3(C1, reg);
            
             A_hat_all_2(1,i,:,:) = our_recovery_adjacency(C2, reg);
             A_hat_all_2(2,i,:,:) = our_recovery_adjacency_v2(C2, reg);
             A_hat_all_2(3,i,:,:) = our_recovery_adjacency_v3(C2, reg);
%             
%             A_hat_all_3(1,i,:,:) = our_recovery_adjacency(C3, reg);
%             A_hat_all_3(2,i,:,:) = our_recovery_adjacency_v2(C3, reg);
%             A_hat_all_3(3,i,:,:) = our_recovery_adjacency_v3(C3, reg);
            %A_hat_all(4,i,:,:) = our_recovery_adjacency_v4(C, reg);

            i = i+1;
            toc
        end
    end
end

%
A_hat_all = A_hat_all_2;
sg = 7:48;
%sg = 7:48;
rng = [1,2,3];
load('st_data\ave_real_data_prediction.mat')
clrs = {'b','g','k','r'};

mi_mat = mi_matrices{nG};
di_mat = di_matrices{nG};
ndmi_mat = mi_nd_matrices{nG}; 
nddi_mat = di_nd_matrices{nG};
Other_models{1} = mi_mat(sg,sg);
Other_models{2} = di_mat(sg,sg);
Other_models{3} = ndmi_mat(sg,sg);
Other_models{4} = nddi_mat(sg,sg);

for t = 1:numel(epsilons)
    for b = rng
        A_hat = squeeze(A_hat_all(b,t,:,:));
        %****************
        % eliminating local trivial predictions
        nn=size(A_hat,1);    
        trunc_mat=zeros(nn,nn);

        for kk=0:4
            if kk~=0
                trunc_mat=trunc_mat+diag(ones(1,nn-kk),kk)+diag(ones(1,nn-kk),-kk);
            else
                trunc_mat=trunc_mat+diag(ones(1,nn-kk),kk);
            end
        end

        trunc_mat=1-trunc_mat;
        contact = A0(sg,sg);
        A_hat_trunc=A_hat.*trunc_mat;
        contact=contact.*trunc_mat;
        for n =1:numel(Other_models)
           Other_models{n} = Other_models{n}.*trunc_mat; 
        end

        for kk=1:nn-8
            if contact(kk,kk+4)==1 & contact(kk+4,kk+8)==1
                contact(kk,kk+8)=0;
            end
        end

        for kk=1:nn-4
            contact(kk,kk+4)=0;
        end

        [res_our,metric1_our,metric2_our,metric3_our]=eval_protein_contact(A_hat_trunc,contact);
        res_our_tot=[];
        res_our_tot=[res_our_tot,res_our'];
        my_ave(t,b,:)=mean(res_our_tot,2);
        res_other_tot = zeros(numel(Other_models),numel(res_our));
        for m =1:numel(Other_models)
            [res_other(m,:),~,~,~]=eval_protein_contact(Other_models{m},contact);
        end
    end
    figure
    plot(res_other(1,:),'b--')

    hold on
    plot(res_other(2,:),'r--')

    hold on
    plot(res_other(3,:),'b')

    hold on
    plot(res_other(4,:),'r')

    %hold on
    %plot(our_ave,'g')
    for b = rng
        hold on 
        plot(squeeze(my_ave(t,b,:)),clrs{b},'LineWidth',1.5)
        xlabel('top predictions')
        ylabel('average fraction of discovered contacts')
        legend('MI','DI','ND','DI+ND','Rec Adj Eig','Rec Adj C','Rec Adj C Hidd')
        title(['Number of nodes: ' num2str(size(A_hat,1)) ' 7:48'])
    end

end

figure();
plot(my_ave(:,:,end))

