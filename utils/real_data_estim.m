function [L_hat,cond,P] = real_data_estim(X,C,model,h,reg)
    if h == 0
        H = 2;
    else
        H = h;
    end
    cond = zeros(1,4);
    P = zeros(20);
    switch model
        case 'GL-SigRep'
            [L_hat,~,~] = GL_SigRep(X,reg);
        case 'GSm-LR'
            [L_hat, ~, ~] = GSm_LR(C,H,reg,false);
        case 'GSm-GL'
            [L_hat, ~, ~] = GSm_GL(C,H,reg,false);
        case 'GSm-St-GL'
            [L_hat,~] = GSm_St_GL(C,reg);
        case 'GSm-St-LR'
            [L_hat,~] = GSm_St_LR(C,reg);
        otherwise
            disp('Error, unknown model')
    end
end