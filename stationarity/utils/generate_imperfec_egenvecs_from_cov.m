function [V_est,Lambda_est,C_est] = generate_imperfec_egenvecs_from_cov(S,R,C_true)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
%S: true shift (N x N)
%R: number of realizations (scalar)
%If C_true==0, then we generate a synthetic covariance inside this funcion If not, we use the matrix C_true
%OUTPUTS: V_est (N x N) critical and mandatory, rest are optional
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N=length(S(1,:)); 

if R<N 
    print('Problem, not enough realizations, sample covariance will be ill-conditioned')
end

if C_true == 0
    h=rand(1,5);
    H=h(1)*eye(N)+h(2)*S+h(3)*S^2+h(4)*S^3+h(5)*S^4;
    C_true = H*H';
end

X_iid=randn(N,R);
X = sqrtm(C_true)*X_iid;
C_est = (1/R)*(X*X');
[V_est, Lambda_est] = eig(C_est);

