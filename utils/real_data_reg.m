function reg = real_data_reg(model)

    switch model
        case 'GL-SigRep'
            reg.alpha = 0.0032;%1e-3;
            reg.beta = 0.1;
            reg.th = 1e-4; %1e-2;
        case 'GSm-St-GL'
            %norm C
            reg.alpha = 1;
            reg.beta = 0.0081;
            reg.lambda = 0.003;
            reg.gamma1 = 10;
            reg.gamma2 = 1;
            reg.epsilon = 0.9;
        case 'GSm-St-LR' 
            reg.alpha =1;
            reg.beta = 0.0081;
            reg.lambda = 8.11e-6;
            reg.gamma1 = 0.12;
            reg.epsilon = 0.9;
        case 'GSm-GL'
            reg.gamma1 = 1e-3;
            reg.beta = 1e-6;
            reg.lambda = 0.1;
            reg.gamma2 = 1;
        case 'GSm-LR'
            reg.gamma1 = 1;
            reg.beta = 4e-3;
            reg.lambda = 1.3e-3;
        otherwise
            disp('Error, unknown model')
    end
end