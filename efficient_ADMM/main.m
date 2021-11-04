addpath('..\utils');
N = 20;
h = 2;
O = N-h;
prms.N = N;
prms.norm_L = false;
prms.L_bin = true;
prms.p = 0.2;
prms.M = 100;
prms.sigma = 0.1;
prms.norm_noise = true;
prms.sampled = true;
g_type = 'ER';
sig_type = 'FA';
links_type = 'rand';
verbose = false;

%Generate graph
[A,L] = generate_graph(g_type,prms);
%L = L/trace(L)*N;
%Generate signals
[X,~,C,snr] = generate_graph_signals(sig_type, L, prms, verbose);
%Select hidden nodes
[s_n, s_h] = select_hidden_nodes(links_type, N-h, L, X);


%Smoothness with A
A_aux = zeros(N,N);
D = zeros(N,N);
for i = 1:N
    for j = 1:N
        D(i,j) =  norm(X(i,:)-X(j,:))^2;
        A_aux(i,j) = A(i,j)*D(i,j);
    end
end
A_smth = sum(sum(A_aux));

%Smoothness with A hadamard prouduct
A_smth_had = norm(vec(A.*D),1);

%Smoothness with L
L_smth = 2*trace(X'*L*X);


%Hidden variables
Ao = A(s_n,s_n);
Aoh = A(s_n,s_h);
Aho = Aoh';
Ah = A(s_h,s_h);

Do = D(s_n,s_n);
Doh = D(s_n,s_h);
Dh = D(s_h,s_h);

Xo = X(s_n,:);
Xh = X(s_h,:);

%norma 1 de A*D
A_smth_had_hid = norm(vec(Ao.*Do),1) + 2*norm(vec(Aoh.*Doh),1) + norm(vec(Ah.*Dh),1);

% Logaritmo de A
A_sum = sum(log(sum(A)));
A_sum_hidd = sum(log(sum(Ao,2)+ sum(Aoh,2))) + sum(log(sum(Aho,2)+ sum(Ah,2)));

%Do, Doh y Dh en funcion de Xo y Xh 
Dox = zeros(O,O);
for i = 1:O
    for j = 1:O
        Dox(i,j) =  norm(Xo(i,:)-Xo(j,:))^2;
    end
end

figure()
subplot(221)
imagesc(Do)
colorbar()
subplot(222)
imagesc(Dox)
colorbar()
subplot(223)
imagesc(Do-Dox)
colorbar()
subplot(224)
imagesc(Do)
colorbar()
%% Version vectorizada
At = A.';
m  = triu(true(size(At)));
a  = At(not(m)).';

Dt = D.';
m  = triu(true(size(Dt)));
d  = Dt(not(m)).';

A_smth_vec = 2*d*a';


% Tensor para calcular la smooothness 
%X1 = repmat(X,[1,1,20]);
%X2 = repmat(X',[1,1,20]);
%X2 = permute(X2,[2,1,3]);





