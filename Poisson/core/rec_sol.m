function out = rec_sol(gf_mat, BF, omega, C, K)

ind = 1:K;
out = 0*omega;

for k = ind
    out = out + C(k)*squeeze(BF(k,:,:)).*omega;
end % for k

out = out + gf_mat;