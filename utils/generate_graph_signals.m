function [X,X_0,C,snr] = generate_graph_signals(sig_type, L, params, verbose)
    if nargin < 4
        verbose = false;
    end
    N = size(L,1);
    M = params.M;
    sigma = params.sigma;
    h_mu = zeros(N,1);
    S = diag(diag(L)) - L;
    switch sig_type
        case 'FA' % Factor analysis
            [V, Lambda] = eig(L); 
            Lambda_inv = pinv(Lambda);
            H = mvnrnd(h_mu, Lambda_inv, M)';
            X_0 = V*H;
            C = pinv(L);% + sigma*I;
        case 'FoE2' % Function of square eigenvalues
            [V, Lambda] = eig(L);
            size(Lambda);
            Sigma_frecs = diag(1./(0.5 + diag(Lambda).^2));
            H = mvnrnd(h_mu, Sigma_frecs, M)';
            X_0 = V*H;
            C = V*Sigma_frecs*V';
        case 'FoE' % Function of eigenvalues
            [V, Lambda] = eig(L);
            size(Lambda);
            Sigma_frecs = diag(1./(0.5 + diag(Lambda)));
            H = mvnrnd(h_mu, Sigma_frecs, M)';
            X_0 = V*H;
            C = V*Sigma_frecs*V';   
        case 'ST'    
            h1 = randn(5,1); % Draw the coefficients of the first polynomial
            H1 = zeros(N,N);
            for ii = 1:5
                H1 = H1 + h1(ii)*S^(ii-1);
            end
            C = H1^2;
            X_iid=randn(N,M);
            X_0 = sqrtm(C)*X_iid;
        otherwise
            error('ERR: Unkown signal type')
    end
    
    p_x = norm(X_0, 'fro')^2/M;
    
    % Set sigma for normalizing the noise power
    if isfield(params, 'norm_noise') && params.norm_noise
        sigma = sqrt(sigma*p_x/N);
    end
    noise = randn(N, M)*sigma;
    
    X = X_0 + noise;
    if params.sampled
       C = X*X'/M;
    end
    p_n = norm(noise, 'fro')^2/M;
    snr = p_x/p_n;
    if verbose
        eq_sigma = sigma^2*N/p_x;
        disp(['Mean Np: ' num2str(p_n) '   norm sigma: ' num2str(eq_sigma)])
        disp(['SNR(nat|dB): ', num2str(snr) ' | ' num2str(10*log10(snr))])
        disp(['Mean smoothness: ' num2str(trace(X_0'*L*X_0)/M)])
    end
end
