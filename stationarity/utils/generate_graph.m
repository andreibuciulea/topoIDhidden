function [A, L] = generate_graph(g_type,params)   
    N = params.N;
    norm_L = params.norm_L;
    sigma = params.sigma;
    max_tries = 10;
    A = [];
    L = [];
    switch g_type
        case 'ER'
            p = params.p;
            A = generate_connected_ER(N,p);
        case 'TREE'
            A = preferential_attachment_graph(N,1);
        case 'BA'
            m = params.m;
            %A = scale_free(N,m,m0);
            A = preferential_attachment_graph(N,m);
        case 'RBF'
            connected = params.connected;
            T = 0.75;
            s = 0.5;
            A = gaussian_graph(N,T,s);
            if connected %generate connected graph
                k = 1;
                lambdas = eig(diag(sum(A))-A);
                while((abs(lambdas(2)) < 1e-6) && (k < max_tries))
                    A = gaussian_graph(N,T,s);
                    k = k+1;
                    lambdas = eig(diag(sum(A))-A);
                end
                if k >= max_tries
                    disp('ERR: Generated unconnected RBF graph')
                end
            end
        case 'RING'
            K = params.K;
            A = small_world(N,K,0);
        case 'SW'
            K = params.K;
            beta = params.beta;
            A = small_world(N,K,beta);
%             mean(sum(A))
        case 'SBM'
            disp('NOT AVAILABLE YET')
        otherwise
            disp('ERR: Unkown graph type')
            return
    end
    
    L = diag(sum(A,2)) - A;
    if norm_L
        %L = L/trace(L)*N;
        L = L/norm(L,'fro');
    else
        sigma = sigma*sqrt(N/trace(L));
    end