function [fsc,out] = get_results(all_prms)
    N = all_prms.N;
    nG = all_prms.nG;
    H = all_prms.H;
    g_type = all_prms.g_type;
    sig_type = all_prms.sig_type;
    rm_node = all_prms.rm_node;
    models = all_prms.models;
    iters = all_prms.iters;
    verb = all_prms.verb;
    fsc = zeros(nG,numel(rm_node),numel(H),numel(models));
    all_norms = zeros(nG,numel(rm_node),numel(H),numel(models));
    smooth = zeros(nG,1);
    parfor g = 1:nG
        %display
        %disp(['Graph: ' num2str(g)])
        %generate graph
        [A, L] = generate_graph(g_type,all_prms);
        %generate signals
        [X,~,C] = generate_graph_signals(sig_type,L,all_prms,verb);
              
        fsc_g = zeros(numel(rm_node),numel(H),numel(models));
        norm_g = zeros(numel(rm_node),numel(H),numel(models));
        for k = 1:numel(rm_node)
            %node selection
            links_type = rm_node{k};
            %display
            %disp(['  Node selection: ' links_type])
            n = 1;
            for h = H
                %display
                %disp(['    Hidden = ' num2str(h)])
                %select hidden nodes
                [s_n, s_h] = select_hidden_nodes(links_type, N-h, L, X);
                Xo = X(s_n,:); Ao = A(s_n,s_n); Co = C(s_n,s_n); Lo = L(s_n,s_n);
                for m = 1:numel(models)
                    %Get fscore for all models
                    model = models{m};
                    reg = get_reg(g_type,sig_type,model);
                    %disp(['      Model: ', model])
                    reg.S0 = Ao;
                    [L_hat, ~] = estimate_graph(model, Xo, Co, h, reg, iters, verb);
                    A_aux = diag(diag(L_hat))-L_hat;
                    Aom = Ao/max(max(Ao));
                    norm_g(k,n,m) = norm(Aom-A_aux/max(max(A_aux))^2,'fro')/norm(double(Aom)^2,'fro');
                    A_hat = mbinarize(A_aux,2);
                    Ao = mbinarize(Ao,2);
                    [~,~,fsc_g(k,n,m),~,~] = graph_learning_perf_eval(Ao,A_hat);
                    
                end
                n = n+1;
            end
        end
        smooth(g) = vec(Co)'*vec(Lo)/norm(Lo);
        fsc(g,:,:,:) = fsc_g;
        all_norms(g,:,:,:) = norm_g;
    end
out.smooth = smooth;
out.all_norms = all_norms;
end