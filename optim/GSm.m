function [L] = GSm(C,reg)
    N = size(C,1);
    alpha = reg.alpha;
    beta = reg.beta;
    lambda = reg.lambda;
    
    %disp(['alpha:' num2str(alpha) '  beta:' num2str(alpha) '  lambda:' num2str(alpha)])
    cvx_begin quiet
    variable L(N,N) symmetric semidefinite
    minimize(alpha*vec(C)'*vec(L) + beta*vec(L(~eye(N)))'*vec(L(~eye(N)))...
        - lambda*ones(1,N)*log(diag(L)))
    subject to
        sum(abs(L*ones(N,1))) <= 1e-8;
        L(~eye(N)) <= 0;
        %L(~eye(N)) >= -1;
    cvx_end
end










