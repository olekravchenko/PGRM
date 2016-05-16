function output = ua(id,x,y)

switch id
    case '01'
        output = sin(pi*x).*cos(2*pi*y);
    case '02'
        s = 0;
        for k = 1:100
            kk = 2*k - 1;
            s = s + ((sinh(kk*pi*(1-y)) + sinh(kk*pi*y)) / sinh(kk*pi)) .* (sin(kk*pi*x) / kk^3);
        end
        output = 0.5 * (x-1) .* x + (4/pi^3) * s;        
end