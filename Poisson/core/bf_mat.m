function out = bf_mat(ind_mat, K, nx, ny, x, y, hf, nf)

out     = zeros(K,ny,nx);
ind     = 1:K;

% compute matrix of basis functions
for k = ind
    i = ind_mat(2,k);
    j = ind_mat(3,k);
    tmp_01 = bdspln_function(i-1, (x - x(1,1))/(x(nx,nx) - x(1,1)), hf, nf);
    tmp_02 = bdspln_function(j-1, (y - y(1,1))/(y(nx,nx) - y(1,1)), hf, nf);
    tmp_03 = tmp_01.*tmp_02;
    out(k,:,:)   = tmp_03;
end