function [bc_mat, rhs_mat, gf_mat] = init_mat(id, x, y)



switch id
    case 1
        bc_mat  = 0*x;
        rhs_mat = 0*y - 1;         
        gf_mat  = bc_mat;
    case 2
        bc_mat      = 0*x;        
        rhs_mat     = 0*y - 1; 
        bc_mat(1,:) = 1;
        bc_mat(:,1) = 1;
        
        tmp1                = r_con(x, y);
        tmp2                = r_con(1 - x, 1 - y);
        gf_mat              = tmp2 ./ (tmp1 + tmp2);
%         gf_mat(gf_mat<0)    = 0;
    case 3
        bc_mat  = 0*x;        
        rhs_mat = -8*pi^2*sin(2*pi*x)*sin(2*pi*y);         
        gf_mat  = bc_mat;
    otherwise
        bc_mat  = 0*x;        
        rhs_mat = 0*y - 1; 
        gf_mat  = bc_mat;
end % switch id
        
