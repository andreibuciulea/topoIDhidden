function [S_hat,out] = GSR_H_Rw_Fact(C,reg,verbose)

N = size(C,1);
tau = reg.tau;
del = reg.del;
max_iter = reg.max_iter;
alpha = reg.alpha;
beta = reg.beta;
H = reg.H;
if verbose
    disp(' -Initialization: Low Rank optimization');
end
reg0 = get_reg(reg.g_type,reg.Ctype,'GSt-Rw',reg.max_iter);
[S_init,out0] = GSR_H_Rw(C,reg0,verbose);

[U,Sigma,V]=svd(out0.P);
F_hat = U(:,1:H)*Sigma(1:H,1:H);
G_hat = V(:,1:H)';

G_prev = G_hat;
F_prev = F_hat;

if isnan(S_init)
    S_prev = ones(N);
else
    S_prev = S_init;
end

if verbose
    disp('  Starting Low Rank reweighted optimization...');
end
for i=1:max_iter
    % Estimate S^(i) from F^(i-1) and G^(i-1)
    weigh_S = tau*ones(N)./(abs(S_prev) + del*ones(N));
    cvx_begin quiet
        variable S_hat(N,N) symmetric
        minimize (alpha*norm(C*S_hat+F_hat*G_hat-S_hat*C-G_hat'*F_hat','fro') + weigh_S(:)'*S_hat(:))
        subject to
        S_hat >= 0;
        diag(S_hat) <= 10^-6;
        abs(sum(S_hat(:,1)) - 1) <= 1e-6;
        %S_hat*ones(N,1) >= 1;
    cvx_end
   
    
    if isnan(cvx_optval)
        disp('WARNING: S: unfeasible problem')
        S_hat = S_prev;
    else
        S_prev = S_hat;
    end
    
    % Estimate G^(i) from F^(i-1) and S^(i)
    weigh_G = tau*ones(H,N)./(abs(G_prev)+ del*ones(H,N));
    cvx_begin quiet
        variable G_hat(H,N)
        minimize (beta*norm(C*S_hat+F_hat*G_hat-S_hat*C-G_hat'*F_hat','fro') + norm(weigh_G(:)'*G_hat(:),1))
        subject to
        G_hat >= 0;
        G_hat <= 1;
    cvx_end
    if isnan(cvx_optval)
        G_hat = G_prev;
    else
        G_prev = G_hat;
    end
    
    % Estimate F^(i) from G^(i) and S^(i)
    cvx_begin quiet
        variable F_hat(N,H)
        minimize (norm(C*S_hat + F_hat*G_hat - S_hat*C - G_hat'*F_hat', 'fro'))

    cvx_end
    if isnan(cvx_optval)
        disp('WARNING: F: unfeasible problem')
        F_hat = F_prev;
    else
        F_prev = F_hat;   
    end 
    %disp(norm(C*S_hat + F_hat*G_hat - S_hat*C - G_hat'*F_hat', 'fro'))
end
out.F_hat = F_hat;
out.G_hat = G_hat;



end