function reg = get_reg(g_type,Ctype,model,max_iter)
    options = [g_type '-' Ctype '-' model];
    
    if strcmp(options,'ER-Cpoly-GSt')
        %reg = struct('eps1',1e-6,'max_iter',max_iter);
        reg = struct('alpha',1e-4,'beta',10,'eps1',1e-6,'max_iter',max_iter);
    elseif strcmp(options,'ER-Cmrf-GSt')
        %reg = struct('eps1',1e-6,'max_iter',max_iter);
        reg = struct('alpha',2.5e-4,'beta',1,'eps1',1e-6,'max_iter',max_iter);
        
    elseif strcmp(options,'ER-Cpoly-GSt-Rw')
        reg = struct('alpha',1e-4,'beta',1,'eps1',1e-6,'tau',1,'del',1e-3,...
            'max_iter',max_iter,'g_type',g_type,'Ctype',Ctype);
    elseif strcmp(options,'ER-Cmrf-GSt-Rw')
        reg = struct('alpha',1.5e-4,'beta',1,'eps1',1e-6,'tau',1,'del',1e-3,...
            'max_iter',max_iter,'g_type',g_type,'Ctype',Ctype);
        
    elseif strcmp(options,'ER-Cpoly-GSt-Rw-Fact') 
        reg = struct('alpha',4.64,'beta',1e3,'tau',1,'del',1e-3,...
            'max_iter',max_iter,'g_type',g_type,'Ctype',Ctype);
    elseif strcmp(options,'ER-Cmrf-GSt-Rw-Fact')
        reg = struct('alpha',1e5,'beta',1e3,'tau',1,'del',1e-3,...
            'max_iter',max_iter,'g_type',g_type,'Ctype',Ctype);
    
        
    elseif strcmp(options,'ER-Cmrf-LVGL') 
        reg = struct('alpha',1e-3,'beta',1e-3);
    elseif strcmp(options,'ER-Cpoly-LVGL') 
        reg = struct('alpha',1,'beta',1e-3);
    else
        reg = {};
        error(['Error: unknown reg for ' options])
    end
    
end