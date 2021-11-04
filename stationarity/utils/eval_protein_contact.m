function [sum_enrich,metric1,metric2,metric3]=eval_protein_contact(net,contact)

% sum_enrich shows the total discovery rate
% metric1 :non-redundant discovery rate (no indirect path)
% metric2: non-redundant iscovery rate, window based (w=1)
% metric3: non-redundant iscovery rate, window based (w=2)

contact=(contact>0)*1.0;
net=triu(net);
contact=triu(contact);
n=size(contact,1);
net_vec=net(1:end);
contact_vec=contact(1:end);

%*********
%sorting predictions
[~,I]=sort(net_vec,'descend');

% only consider top 250 predictions
net_vec_sorted=net_vec(I(1:250));
contact_vec_sorted=contact_vec(I(1:250));
l=length(net_vec_sorted);
sum_enrich=zeros(1,l);

metric1=zeros(1,l); % non-redundant
metric2=zeros(1,l); % w=1
metric3=zeros(1,l); % w=2


contact_temp1=contact*0;
contact_temp2=contact*0;
contact_temp3=contact*0;

% Initialization of metrics for the first prediction
if contact_vec_sorted(1)==1
    sum_enrich(1)=1;
    [x,y]=index_to_pair(I(1),n);
    
    metric1(1)=1;
    metric2(1)=1;
    metric3(1)=1;
    
    contact_temp1(x,y)=1;
    contact_temp1(y,x)=1;
    
    contact_temp2(x,y)=1;
    contact_temp2(y,x)=1;
    
    contact_temp3(x,y)=1;
    contact_temp3(y,x)=1;
end

for i=2:l
    %****
    % total discovery rate
    sum_enrich(i)=sum_enrich(i-1)+contact_vec_sorted(i);
    
    
    [x,y]=index_to_pair(I(i),n);
    %****
    % non-redundant (no path)
    b=graphshortestpath(sparse(contact_temp1),x,y);
    if b==inf & contact_vec_sorted(i)==1
        
        metric1(i)=metric1(i-1)+1;
    else
        metric1(i)=metric1(i-1);
    end
    
    contact_temp1(x,y)=1;
    contact_temp1(y,x)=1;
    %*********
    % non-redundant (window based)
    
    w=1;
    x_low=max(x-w,1);
    x_high=min(x+w,n);
    y_low=max(y-w,1);
    y_high=min(y+w,n);
    
    if sum(sum(contact_temp2(x_low:x_high,y_low:y_high)))==0 & contact_vec_sorted(i)==1
        metric2(i)=metric2(i-1)+1;
    else
        metric2(i)=metric2(i-1);
    end
    
    if contact_vec_sorted(i)==1
        contact_temp2(x,y)=1;
        contact_temp2(y,x)=1;
    end
    
    %*********
    % non-redundant (window based)
    
    w=2;
    x_low=max(x-w,1);
    x_high=min(x+w,n);
    y_low=max(y-w,1);
    y_high=min(y+w,n);
    
    if sum(sum(contact_temp3(x_low:x_high,y_low:y_high)))==0 & contact_vec_sorted(i)==1
        metric3(i)=metric3(i-1)+1;
    else
        metric3(i)=metric3(i-1);
    end
    
    if contact_vec_sorted(i)==1
        contact_temp3(x,y)=1;
        contact_temp3(y,x)=1;
    end
 
end

sum_enrich=sum_enrich./sum(sum(contact));
metric1=metric1./sum(sum(contact));
metric2=metric2./sum(sum(contact));
metric3=metric3./sum(sum(contact));







