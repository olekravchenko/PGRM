function [out1, out2] = GetFunc(id)

if id == 1
    rhs_func    = @(x,y) 2*exp(x+y).*(-1-2*y+x.^2.*y.*(2+y)+2*x.*(y.^2-1));    % Right hand side function
    ua          = @(x,y) (1-x.^2).*(1-y.^2).*exp(x+y);      % Exact Solution
elseif id == 2
    rhs_func    = @(x,y) -2*pi*sin(x).*sin(y);    % Right hand side function
    ua          = @(x,y) sin(pi*x).*sin(pi*y);      % Exact Solution
elseif id == 3
    ua          = @(x,y) (1 - x.^2).*(1 - y.^2);
    rhs_func    = @(x,y) 2*(x.^2 + y.^2 - 2);
else
    disp('Incorrect Input');
    out1 = 0; out2 = 0;
    return
end % if

out1 = ua;
out2 = rhs_func;