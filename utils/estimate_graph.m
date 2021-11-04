function [L_hat, R_hat] = estimate_graph(model, Xo, Co, H, reg, iters, verb)
    if nargin < 6
        verb = false;
    end
    R_hat = [];
    S_hat = 0;
    reg.H = H;
    switch model
        case 'GL-SigRep'
            [L_hat,~,~] = GL_SigRep(Xo,reg);
        case 'GSm'
            [L_hat] = GSm(Co,reg);  % Xo = Co
        case 'GSm-LR'
            [L_hat, ~, R_hat] = GSm_LR(Co,H,reg,verb);
        case 'GSm-GL'
            [L_hat, ~, R_hat] = GSm_GL(Co,H,reg,verb);
        case 'GSm-St-LR'
            [L_hat, ~] = GSm_St_LR(Co,reg);
        case 'GSm-St-GL'
            [L_hat, ~] = GSm_St_GL(Co,reg);
        case 'GSt'
            %[L_hat, ~, R_hat] = GSt(Co,H,reg,verb);
        case 'GSt-Rw-Fact'
            %[L_hat, ~, R_hat] = GSt_Rw_Fact(Co,H,reg,verb);
        case 'LVGL'
            %[L_hat, ~, R_hat] = LVGL(Co,H,reg,verb);      
        otherwise
            error('ERR: Unkown estimation method')
    end
    if sum(sum(S_hat)) ~= 0
        L_hat = diag(diag(S_hat)) - S_hat; 
    end
end
