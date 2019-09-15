
%% To pass examples in each tree for getting root-leaf paths, from the root
%

function [rl_paths, path_lengths] = channel_one_tree(X, classRF_model, idx_tree, max_path_len)   


    NODE_TERMINAL=-1;
%     NODE_TOSPLIT =-2;
%     NODE_INTERIOR=-3;
    N = size(X,1);
    
    % paths from root to a leaf
    rl_paths = zeros(N, max_path_len);
    path_lengths = zeros(N, 1);
    for idx_samp = 1 : N
        k = 1;
        path_ins = ones(1, max_path_len) * -10 * idx_samp;
        path_ins(1) = 1;
        % start with k=1 and then go on till we reach a terminal node
        % nodestatus is i think numnodes x numberof trees
        path_len = 1;
        while (classRF_model.nodestatus(k,idx_tree) ~= NODE_TERMINAL)
    	    % m is the variable that was used to split
            m = classRF_model.bestvar(k,idx_tree);
            
            % now that we know m we can find if X(current_example,m) <= the split value
            %then either go right or left in the tree.           
            % if the X's value is less then the split go left, else go right
            if X(idx_samp,m)  <= classRF_model.xbestsplit(k,idx_tree)
                %k = treemap((k-1)*2+1);
                k = classRF_model.treemap2((k-1)*2+1,idx_tree);
                
            else
                %k = treemap((k-1)*2+2);
                k = classRF_model.treemap2((k-1)*2+2,idx_tree);
                
            end
            
            path_len = path_len + 1;
            path_ins(path_len) = k;
            
        end
        
        path_lengths(idx_samp) = path_len;
        rl_paths(idx_samp, :) = path_ins;   
            
    end
end


