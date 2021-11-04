function [Lo, K, r] = GSm_LR(Co,H,reg,verbose)
    O = size(Co,1);
    gamma1 = reg.gamma1;
    beta = reg.beta; 
    lambda = reg.lambda;
    if nargin < 4
        verbose = false;
    end
    cvx_begin quiet
        variable Lo(O,O) symmetric semidefinite
        variable K(O,O) 
        variable r
        minimize(vec(Co)'*vec(Lo)+2*trace(K)+ r ...
            + beta*Lo(~eye(O))'*Lo(~eye(O))...
            - lambda*ones(1,O)*log(diag(Lo))...
            + gamma1*norm_nuc(K))
        subject to
            vec(Co)'*vec(Lo)+2*trace(K)+r >= 0;
            r >= 0;
            Lo(~eye(O)) <= 0;
            sum(abs(sum(Lo))) <= 1e-8;
    cvx_end
    if verbose
        disp(['EST: tr(CoLo): ' num2str(vec(Co)'*vec(Lo))...
            '  tr(K): ' num2str(trace(K)) '  r: ' num2str(r)])
    end
end