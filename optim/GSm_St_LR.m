function [L,out] = GSm_St_LR(C,reg)
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
        variable L(N,N) symmetric semidefinite
        variable L_hat(N,N) symmetric semidefinite
        variable Q(N,N) 
        variable R(H,H) symmetric semidefinite
        minimize(alpha*(vec(C)'*vec(L)+2*trace(Q)+trace(R))...
            + beta*L(~eye(N))'*L(~eye(N))...
            - lambda*ones(1,N)*log(diag(L))...
            + gamma1*norm_nuc(Q))
        subject to
            vec(C)'*vec(L)+2*trace(Q)+trace(R) >= 0;
            L(~eye(N)) <= 0;
            %sum(L) == 0;
            sum(abs(sum(L))) <= 1e-8;
            norm(C*L_hat-L_hat*C+Q-Q','fro') <= 1e-6;
            norm(L-L_hat,'fro') <= epsilon;
    cvx_end
    out = 0;
end