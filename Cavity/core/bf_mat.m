function out = bf_mat(ind_mat, K, nx, ny, x, y, hf, nf)

out     = zeros(K,nx,ny);
ind     = 1:K;

% compute matrix of basis functions
for k = ind
    i = ind_mat(2,k);
    j = ind_mat(3,k);
    out(k,:,:)   = bdspln_function(i-1, x, hf, nf).*bdspln_function(j-1, y, hf, nf);
end