function [s_n, s_h] = select_hidden_nodes(sel_type, O, L, C)
    N = size(L,1);
    s_n = [];
    s_h = [];
    A = abs(L);
    A(~~eye(N)) = 0;
    switch sel_type
        % select nodes with less links
        case 'min' 
            connections = sum(A);
            [~,node_pos] = sort(connections,'descend');
        % select nodes with more links
        case 'max'
            connections = sum(A);
            [~,node_pos] = sort(connections,'ascend');
        % select nodes which when removed the smoothness is smaller
        case 'min_smooth'
            smooth = node_smoothness(N, L, C);
            [~,node_pos] = sort(smooth,'descend');
        % select nodes which when removed the smoothness is bigger
        case 'max_smooth'
            smooth = node_smoothness(N, L, C);
            [~,node_pos] = sort(smooth,'ascend');
        % select nodes which when removed the variation of the smoothness
        % is bigger
        case 'max_diff_sm'
            smooth = node_smoothness(N, L, C);
            diff_sm = abs(smooth-trace(C*L));
            [~,node_pos] = sort(diff_sm,'ascend');
        case 'rand'
            [~,node_pos] = sort(rand(N,1));
        % select nodes at random
        otherwise
            disp('ERR: unknown method for selecting the hidden node')
            return
    end
    s_n = sort(node_pos(1:O));
    s_h = sort(node_pos(O+1:N));
end


function smoothness = node_smoothness(N, L, C)
    smoothness = zeros(N,1);
    for i=1:N
        Lo = L([1:i-1 i+1:N], [1:i-1 i+1:N]);
        Co = C([1:i-1 i+1:N], [1:i-1 i+1:N]);
        smoothness(i) = trace(Co*Lo);
    end
end