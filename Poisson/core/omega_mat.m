function omega = omega_mat(id, x, y)

switch id
    case 1 % 5-gone
        om1 = 0.5*(0.5 - (x - 0.5).^2 - (y - 0.5).^2 - ...
            sqrt(x.^2.*(x-1.0).^2 + y.^2.*(y-1.0).^2));
        om2 = 1.5 - x - y;
        om3 = y - 2*x + 1.5;
        om12 = 0.5*(om1 + om2 - sqrt(om1.^2 + om2.^2));
        omega = 0.5*(om12 + om3 - sqrt(om12.^2 + om3.^2));
    case 2 % square
        omega = 0.5*(0.5 - (x - 0.5).^2 - (y - 0.5).^2 - ...
            sqrt(x.^2.*(x-1.0).^2 + y.^2.*(y-1.0).^2));
    case 3 % uniform one        
        omega = 1 + 0*x;
    case 4 % rectangle [a b c d] = [0 b 0 d]
        b = x(end); d = y(end);        
        f1 = x.*(b-x);  f2 = y.*(d-y);
        omega = 0.5*(f1 + f2 - sqrt((f1 - f2).^2));
    otherwise
        omega = 2 - x^2 + y^2; 
end