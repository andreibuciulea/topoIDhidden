function reg = get_reg(g_type,sig_type,model,M)
    options = [g_type '-' sig_type '-' model];
    
    % Reg weights for signal models: RBF-FA
    if strcmp(options,'ER-FA-GL-SigRep')
         reg = struct('alpha',0.0032,'beta',0.1,'max_iter',50,'th',12.59);
         
    elseif strcmp(options,'RBF-FA-GL-SigRep')
         reg = struct('alpha',0.012,'beta',0.79,'max_iter',50,'th',0.06);% L_norm=false, M=100 %exp1
         %reg = struct('alpha',8.86e-4,'beta',0.6,'max_iter',50,'th',0.1);%L_norm=true, M=500
         
    elseif strcmp(options,'RBF-FA-GSm')
        reg = struct('alpha',1,'beta',4.5204,'lambda',0.4175);%exp1
        
    elseif strcmp(options,'RBF-FA-GSm-LR')
        %reg = struct('beta',3.0392,'lambda',0.127,'gamma1',1);%M=100 %exp1
        %reg = struct('gamma',1,'beta',0.0162,'lambda',0.4833);%Perf Cov 
        reg = struct('gamma1',1,'beta',12.9,'lambda',0.0215);%Perf Cov exp7
        
    elseif strcmp(options,'RBF-FA-GSm-GL') 
        %reg = struct('gamma1',1,'beta',0.01334,'lambda',1.333,'gamma2',0); 
        reg = struct('beta',1,'lambda',1,'gamma1',3);%exp7;
        %reg = struct('beta',1.72,'lambda',1,'gamma1',141); %exp1
        
    elseif strcmp(options,'ER-FA-GSm')
        reg = struct('alpha',1,'beta',0.4833,'lambda',0.0886);
        
    elseif strcmp(options,'BA-FA-GSm')
        reg = struct('alpha',1,'beta',0.4833,'lambda',0.0886);
        
    elseif strcmp(options,'SW-FA-GSm')
        reg = struct('alpha',1,'beta',0.04,'lambda',0.25);
        
    
        
    elseif strcmp(options,'ER-FA-GSm-LR')
        %reg = struct('beta',0.1,'lambda',0.1778,'gamma1',1);
        reg = struct('beta',1,'lambda',0.215,'gamma1',1);%exp2
        
    elseif strcmp(options,'ER-FA-GSm-GL')
        %reg = struct('beta',0.01334,'lambda',1.333,'gamma1',1);
        reg = struct('beta',1.29,'lambda',1,'gamma1',10);%exp2 
        
    elseif strcmp(options,'RBF-FA1-GSm-LR') 
        %reg = struct('gamma1',1,'beta',0.01334,'lambda',1.333,'gamma2',0);
        reg = struct('beta',1,'lambda',1e-2,'gamma1',1); %exp3
        
    elseif strcmp(options,'BA-FA-GSm-LR')
        reg = struct('beta',0.1,'lambda',0.1778,'gamma1',1);
        
        
    elseif strcmp(options,'RBF-FA-GSm-St-LR')
        %reg = struct('alpha',1,'beta',0.215,'lambda',2.5,'gamma1',8.11,'epsilon',1e-1);% beta= 1e-1,1e0
        reg = struct('alpha',1,'beta',0.215,'lambda',1,'gamma1',3,'gamma2',0,'epsilon',1e-1);
    elseif strcmp(options,'RBF-FA-GSm-St-GL') 
        reg = struct('alpha',1,'beta',0.2,'lambda',2.5,'gamma1',5,'gamma2',0,'epsilon',1e-1);
        
    else
        error(['ERROR: unknown regularizer weigths for ' options])
    end

end