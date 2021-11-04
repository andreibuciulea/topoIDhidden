function [Lo,out] = GSm_St_GL(Co,reg)
    O = size(Co,1);
    H = reg.H;
    alpha = reg.alpha;
    beta = reg.beta;
    lambda = reg.lambda;
    gamma1 = reg.gamma1;
    epsilon = reg.epsilon;
    disp(['alpha:' num2str(alpha) ' beta:' num2str(beta) ' lambda:' num2str(lambda) ' gamma:' num2str(gamma1) ' Epsilon: ', num2str(epsilon)])
    %disp(['           Epsilon: ', num2str(epsilon)])
    
    cvx_begin quiet
        variable Lo(O,O) symmetric semidefinite
        variable Lo_hat(O,O) symmetric semidefinite
        variable K(O,O)
        variable r
        minimize(alpha*vec(Co)'*vec(Lo)+2*trace(K)+r...
            + beta*Lo(~eye(O))'*Lo(~eye(O))...
            - lambda*ones(1,O)*log(diag(Lo))...
            + gamma1*sum(norms(K,2)))
        subject to
            vec(Co)'*vec(Lo)+2*trace(K)+ r >= 0;
            Lo(~eye(O)) <= 0;
            r >= 0;
            sum(abs(sum(Lo))) <= 1e-8;
            norm(Co*Lo_hat-Lo_hat*Co+K-K','fro') <= 1e-6;
            norm(Lo-Lo_hat,'fro') <= epsilon;
    cvx_end
out.K = K;

end