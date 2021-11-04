function [C] = generate_C(S,prms)
    N = prms.N;
    L = prms.Taps;
    type = prms.Ctype;
    R = prms.Samples;
    
    
    if strcmp(type,'Cpoly')
        h1 = randn(L,1); % Draw the coefficients of the first polynomial
        H1 = zeros(N,N);
        for ii = 1:L
            H1 = H1 + h1(ii)*S^(ii-1);
        end
        C = H1^2;
    elseif strcmp(type,'Cmrf')
        [~,D] = eig(S);
        C_inv = (0.01 - min(diag(D)))*eye(N,N) + (0.9 + 0.1*rand(1,1))*S;
        C = inv(C_inv);
        [~,Dc] = eig(C);
        if min(diag(Dc)) < 0
            disp('no es def pos')
            [C,~] = generate_C(S,prms);
        end
    else
        disp('WARNING: unknown covariance type')
        C = [];
        return
    end
    if prms.noise
        [~,~,C] = generate_imperfec_egenvecs_from_cov(S,R,C);
    end
end