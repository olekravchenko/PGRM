function out = GetOmega(id, a, b)

if id == 1
    omega = @(x,y) (a.^2 - x.^2) .* (b.^2 - y.^2);    
else
    disp('Incorrect Input');
    out = nan;
    return
end % if

out = omega;