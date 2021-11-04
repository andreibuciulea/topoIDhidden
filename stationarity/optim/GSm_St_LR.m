function [Lo,out] = GSm_St_LR(C,reg)
    N = size(C,1);
    H = reg.H;
    alpha = reg.alpha;
    beta = reg.beta;
    lambda = reg.lambda;
    gamma1 = reg.gamma1;
    epsilon = reg.epsilon;
    
    %disp(['alpha:' num2str(alpha) ' beta:' num2str(beta) ' lambda:' num2str(lambda) ' gamma:' num2str(gamma)])
    %disp(['           Epsilon: ', num2str(epsilon)])
    cvx_begin quiet
        variable Lo(N,N) symmetric semidefinite
        variable Lo_hat(N,N) symmetric semidefinite
        variable K(N,N) 
        variable r
        minimize(alpha*(vec(C)'*vec(Lo)+2*trace(K)+r)...
            + beta*Lo(~eye(N))'*Lo(~eye(N))...
            - lambda*ones(1,N)*log(diag(Lo))...
            + gamma1*norm_nuc(K))
        subject to
            vec(C)'*vec(Lo)+2*trace(K)+r >= 0;
            Lo(~eye(N)) <= 0;
            r >= 0;
            sum(abs(sum(Lo))) <= 1e-8;
            norm(C*Lo_hat-Lo_hat*C+K-K','fro') <= 1e-6;
            norm(Lo-Lo_hat,'fro') <= epsilon;
    cvx_end
out.K = K;
end