function [Lo,out] = GSm_St_GL(Co,reg)
    O = size(Co,1);
    H = reg.H;
    alpha = reg.alpha;
    beta = reg.beta;
    lambda = reg.lambda;
    gamma1 = reg.gamma1;
    epsilon = reg.epsilon;
    gamma2 = 0;
    cond = zeros(1,4);
    %disp(['alpha:' num2str(alpha) ' beta:' num2str(beta) ' lambda:' num2str(lambda) ' gamma:' num2str(gamma1) ' gamma2:' num2str(gamma2)  ' Epsilon: ', num2str(epsilon)])
    %disp(['           Epsilon: ', num2str(epsilon)])
    
    cvx_begin quiet
        variable Lo(O,O) symmetric semidefinite
        variable L_hat(O,O) symmetric semidefinite
        variable P(O,O)
        variable R(H,H) symmetric semidefinite
        minimize(alpha*(vec(Co)'*vec(Lo)+2*trace(P)+trace(R))...
            + beta*Lo(~eye(O))'*Lo(~eye(O))...
            - lambda*ones(1,O)*log(diag(Lo))...
            + gamma1*sum(norms(P,2)) + gamma2*norm_nuc(P))
        subject to
            vec(Co)'*vec(Lo)+2*trace(P)+trace(R) >= 0;
            Lo(~eye(O)) <= 0;
            %sum(Lo) == 0;
            sum(abs(sum(Lo))) <= 1e-8;
            norm(Co*L_hat-L_hat*Co+P-P','fro') <= 1e-6;
            norm(Lo-L_hat,'fro') <= epsilon;
    cvx_end
out.R = R;
out.cond = cond;

end