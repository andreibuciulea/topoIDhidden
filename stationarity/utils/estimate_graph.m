function [S_hat, out] = estimate_graph(model,Co,reg,verb)
%     if nargin < 6
%         verb = false;
%     end
    switch model
        case 'GSt'
            [S_hat, out]= GSR_H(Co,reg,verb);
        case 'GSt-Rw'
            [S_hat,out] = GSR_H_Rw(Co,reg,verb);
        case 'GSt-Rw-Fact'
            [S_hat, out] = GSR_H_Rw_Fact(Co,reg,verb);
        case 'LVGL'
            [S_hat, out] = LVGLASSO(Co,reg,verb);
        otherwise
            error(['Error: unknown estimation method' model])
    end
end
