function  [res,out] = obtain_est_results(A_hat,contact)
    nn=size(A_hat,1);    
    trunc_mat=zeros(nn,nn);

    for kk=0:4
        if kk~=0
            trunc_mat=trunc_mat+diag(ones(1,nn-kk),kk)+diag(ones(1,nn-kk),-kk);
        else
            trunc_mat=trunc_mat+diag(ones(1,nn-kk),kk);
        end
    end

    trunc_mat=1-trunc_mat;
    %contact = A0(sg,sg);
    A_hat_trunc=A_hat.*trunc_mat;
    contact=contact.*trunc_mat;

    % for n =1:numel(Other_models)
    %    Other_models{n} = Other_models{n}.*trunc_mat; 
    % end

    for kk=1:nn-8
        if contact(kk,kk+4)==1 && contact(kk+4,kk+8)==1
            contact(kk,kk+8)=0;
        end
    end

    for kk=1:nn-4
        contact(kk,kk+4)=0;
    end

    [res,metric1,metric2,metric3]=eval_protein_contact(A_hat_trunc,contact);
    out.m1 = metric1;
    out.m2 = metric2;
    out.m3 = metric3;
    

end