%% Construct the affinity matrix based on the knowledge encoded in the tree hierarchies
%
% @Author: Xiatian (Eddy) Zhu
% @Date: 16 June 2014

function [A] = build_ClustRF_Strct_A(X, RF_model, model_name)

%% Treemap: a variable that holds the tree structure
RF_model.treemap2 = zeros(RF_model.nrnodes*2,RF_model.ntree);
for idx_tree = 1 : RF_model.ntree
    RF_model.treemap2(:,idx_tree) = [RF_model.treemap(:,(idx_tree-1)*2+1); RF_model.treemap(:,(idx_tree-1)*2+2)];
end

%% Iterate through each tree and get their predictions
A = 0; 

for j = 1 : RF_model.ntree

    [tree_paths, path_lengths] = channel_one_tree(X, RF_model, j, 2*size(X, 1));

	inv_nodesize_table = 1 ./ double(RF_model.node_size_table(:, j));
	
    switch model_name
	
		case 'Uniform'
		
			rmv_last_treeNode = 1;
			rmv_treeRoot = 0;
            
        case 'Adaptive'
            
            rmv_last_treeNode = 0;
            rmv_treeRoot = 1;
        
    end
    
    A_tree = calc_path_similarity(tree_paths, path_lengths, inv_nodesize_table, rmv_treeRoot, rmv_last_treeNode, model_name);
    
    A = A + A_tree;

end

%% average
A = A / RF_model.ntree;

%A(eye(size(X, 1)) == 1) = 1;



%% Compute the similarity between paths
function [A_tree] = calc_path_similarity(tree_paths, path_lengths, inv_nodesize_table, rmv_treeRoot, rmv_last_treeNode, model_name)

N = size(tree_paths, 1);

%% remove redundancy dimensions
max_path_len = max(path_lengths);
tree_paths(:, max_path_len+1:end) = [];

P = kron(tree_paths, ones(N, 1));
Q = repmat(tree_paths, N, 1);

P_len = kron(path_lengths, ones(N, 1));
Q_len = repmat(path_lengths, N, 1);
lens = [P_len'; Q_len'];

switch model_name 
    
    case 'Uniform'

        max_len = max(lens);
        basis_path_parts = max_len' - 1; % Do not count the root in

        nodeIdx_on_shared_path_parts = P == Q;
        shared_path_parts = sum(nodeIdx_on_shared_path_parts, 2) - 1; % Do not count the root
    
        A_tree = shared_path_parts ./ basis_path_parts;
        
		clear shared_path_parts basis_path_parts max_len;
 
    case 'Adaptive'  

		nodesize_table_star = [0; inv_nodesize_table];
        
        %% calc the basis path parts 
        % (tree node index) the longer path is used as the basis one
        [~, ind] = max(lens);
        nodeIdx_on_basis_path_parts = P .* repmat((ind == 1)', 1, size(tree_paths, 2));
        nodeIdx_on_basis_path_parts = nodeIdx_on_basis_path_parts + Q .* repmat((ind == 2)', 1, size(tree_paths, 2));
        
        % nodesize on the basis paths
        nodeIdx_on_basis_path_parts(nodeIdx_on_basis_path_parts > 0) = nodeIdx_on_basis_path_parts(nodeIdx_on_basis_path_parts > 0) + 1;
        nodeIdx_on_basis_path_parts(nodeIdx_on_basis_path_parts <= 0) = 1;
        nodesize_on_basis_path_parts = nodesize_table_star(nodeIdx_on_basis_path_parts);
        
        if rmv_treeRoot
           
            nodesize_on_basis_path_parts(:, 1) = 0;
            
        end
        
        basis_path_parts = sum(nodesize_on_basis_path_parts, 2);        
        
        %% calc the shared path parts
        % (node index)
        nodeIdx_on_shared_path_parts = (P == Q) .* P;
        
        if rmv_last_treeNode
            
            % remove the node where samples divert
            nodeIdx_on_shared_path_parts_act = nodeIdx_on_shared_path_parts > 0;
            last_shared_node_pos = sum(nodeIdx_on_shared_path_parts_act, 2);
            last_shared_node_idx = (last_shared_node_pos-1)*size(P, 1) + (1:size(P, 1))';
            nodeIdx_on_shared_path_parts(last_shared_node_idx) = 0;
            
            clear nodeIdx_on_shared_path_parts_act last_shared_node_pos last_shared_node_idx;
            
        end
        
        % nodesize on the shared path parts
        nodeIdx_on_shared_path_parts(nodeIdx_on_shared_path_parts>0) = nodeIdx_on_shared_path_parts(nodeIdx_on_shared_path_parts>0) + 1;
        nodeIdx_on_shared_path_parts(nodeIdx_on_shared_path_parts <= 0) = 1;
        nodesize_on_shared_path_parts = nodesize_table_star(nodeIdx_on_shared_path_parts);
        
        if rmv_treeRoot
           
            nodesize_on_shared_path_parts(:, 1) = 0;
            
        end

        shared_path_parts = sum(nodesize_on_shared_path_parts, 2);

        A_tree = shared_path_parts ./ basis_path_parts;
        
        clear nodesize_on_shared_path_parts nodeIdx_on_shared_path_parts shared_path_parts basis_path_parts;        
    
end

clear P Q P_len Q_len lens;

A_tree = reshape(A_tree, N, N);

switch model_name 

	case 'Uniform'
		% Nothing to do
		
    case 'Adaptive' 
        A_tree = (A_tree + A_tree') / 2;
        
end

%% symmetric or not
if sum(sum(abs(A_tree - A_tree'))) > 0
    error('Not symmetric');
end
    