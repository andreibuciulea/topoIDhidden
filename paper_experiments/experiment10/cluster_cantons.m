function clstr = cluster_cantons(S)
    N = size(S,1);
    S_bin = mbinarize(S,2);%Adj binarizada 
    S_hat = S.*S_bin; % Adj binarizada con pesos
    clstr = zeros(N,2);
    %falta clusterizar los nodos
    L = diag(sum(S_bin))-S_bin;
    [U,~] = eig(L);
    U_aux = U(:,2:4);
    clstr(:,1) = kmeans(U_aux,3);
    %clusterizar por pesos
    L_hat = diag(sum(S_hat))-S_hat;
    [U_hat,~] = eig(L_hat);
    U_aux_hat = U_hat(:,2:4);
    clstr(:,2) = kmeans(U_aux_hat,3);

end