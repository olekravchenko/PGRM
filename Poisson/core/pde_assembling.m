function [A, B] = pde_assembling(init_mat, bf_mat, omega, K, tx, ty, dx, dy)

A   = zeros(K,K);
B   = zeros(1,K);
ind = 1:K;

for l = ind
    B(l)        = trapz(tx, trapz(ty, init_mat.*omega.*squeeze(bf_mat(l,:,:))));    
    for k = ind                
        tmp     = del2(omega.*squeeze(bf_mat(l,:,:)), dx, dy);
        A(k, l) = trapz(tx, trapz(ty, omega.*squeeze(bf_mat(k,:,:)).*tmp));
    end % for l
end % for k