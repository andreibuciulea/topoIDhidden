function reg = sel_reg(model,H)
    epsilon = 1e-2;
    lambda = 1;
    
    if strcmp(model,'GL-SigRep')
        gamma=1;
        alpha = 1e-2;
        if H == 1
            beta = 22;
        elseif H == 2
            beta = 22;
        elseif H == 3
            beta = 16.5; %falta ajusta esto
        end
    elseif strcmp(model,'GSm-St-LR')
        if H == 1
            alpha = 10;beta = 1;gamma = 0.1;
        elseif H == 2
            alpha = 1;beta = 10;gamma = 0.1;
        elseif H == 3
            alpha = 1.5;beta = 10;gamma = 0.1;%falta ajusta esto
        end    
    end
    reg = struct('alpha',alpha,'beta',beta,'lambda',lambda,'gamma1',gamma,...
                        'gamma2',0,'gamma3',0,'epsilon',epsilon,'th',epsilon, 'H',H);
end