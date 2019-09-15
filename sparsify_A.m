%% Make the affinity sparse

function sparse_A = sparsify_A(A)

cur_k = 20;
sparse_A = A;
N = size(sparse_A, 1);
for idx_rwo = 1 : N
    [~, samp_idx] = sort(sparse_A(idx_row, :), 'descend');
    sparse_A(idx_row, samp_idx(cur_k+2:end)) = 0;
end
sparse_A = max(sparse_A, sparse_A'); % make it symmetric
sparse_A(eye(N) == 1) = 0;
