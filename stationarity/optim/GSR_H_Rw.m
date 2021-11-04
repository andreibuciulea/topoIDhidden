function [S_hat,out] = GSR_H_Rw(C,reg,verbose)

N = size(C,1);
alpha = reg.alpha;
beta = reg.beta;
eps1 = reg.eps1;
max_iter = reg.max_iter;
tau = reg.tau;
del = reg.del;

reg0 = get_reg(reg.g_type,reg.Ctype,'GSt',3);
[S_init,out] = GSR_H(C,reg0,verbose);



if isnan(S_init)
    S_prev = ones(N);
else
    S_prev = S_init;
end

for i=1:max_iter
    weigh_S = tau*ones(N)./(abs(S_prev) + del*ones(N));
    cvx_begin quiet
        variable S_hat(N,N) symmetric
        variable P_hat(N,N) 
        minimize (alpha*weigh_S(:)'*S_hat(:) + beta*norm_nuc(P_hat))
        subject to
            abs(diag(S_hat)) <= 1e-6;
            S_hat >= 0;
            S_hat*ones(N,1) >= 1;
            norm(C*S_hat + P_hat - S_hat*C - P_hat','fro') <= eps1;
    cvx_end
   
    %disp(num2str(norm(S_hat-S_prev,'fro')))
    if norm(S_hat-S_prev,'fro') < 1e-3
        %disp('break;')
        break;
    end
    if sum(sum(isnan(S_hat))) ~= 0
        disp('ERROR: unfeasible problem')
        S_hat = S_prev;
    else
        S_prev = S_hat;
    end
    
end
out.scr = 0;
out.P = P_hat;

end