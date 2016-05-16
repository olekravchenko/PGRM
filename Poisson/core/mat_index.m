function ind_mat = mat_index(nf)

ind_mat = zeros(3,(nf+2)^2);

for i = 1:nf+2
    for j = 1:nf+2
        k = (i-1)*(nf+2) + j;
        % fill matrix of indexes
        ind_mat(1,k) = k;
        ind_mat(2,k) = i;
        ind_mat(3,k) = j;        
    end % for j
end % for k