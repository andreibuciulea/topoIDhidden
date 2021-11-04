function [S_hat,out] = GSR_H(C,reg,verbose)
N = size(C,1);
alpha = reg.alpha;
beta = reg.beta;
eps1 = reg.eps1;
max_iter = reg.max_iter;

if verbose
   disp('  -Starting GSR Low Rank optimization...') 
end

for i = 1:max_iter
    cvx_begin quiet
        variable S_hat(N,N) symmetric
        variable K_hat(N,N) 
        %        Sparse matrix  +  low rank matrix
        minimize (alpha*norm(S_hat(:),1) + beta*norm_nuc(K_hat))
        subject to
        diag(S_hat) <= 1e-6;
        S_hat >= 0;
        S_hat*ones(N,1) >= 1;
        norm(C*S_hat + K_hat - S_hat*C - K_hat', 'fro') <= eps1;

    cvx_end
    if (~isnan(sum(sum(S_hat))))
       break;
    else
        disp('CVX optval error')
        eps1 = eps1*2;
    end
end
out.K = K_hat;
if verbose 
    figure(1)
    subplot(2,5,3)
    imagesc(S_hat)
    colorbar()
    title('S GSt')
    subplot(2,5,8)
    imagesc(K_hat)
    colorbar()
    title('K GSt')
end
end






