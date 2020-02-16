function [U_xx U_xy U_xz U_yy U_yz U_zz]  = GetDerivative_6dim(mu, Y) 


    x = Y(1, :);
    y = Y(2, :);
    z = Y(3, :);
    xdot = Y(3, :);
    ydot = Y(4, :);
    zdot = Y(6, :);
    
    d = sqrt((x+mu).^2 + y.^2 + z.^2);
    r = sqrt((x -1 + mu).^2 +y.^2 + z.^2);
        
    U_xx(:) = 1 - (1-mu)./d.^3 - mu./r.^3 + 3*(1-mu)./d.^5.*(x+mu).^2 + 3.*mu./r.^5.*(x-1+mu).^2;
    U_xy(:) = 3*(1-mu)./d.^5.*(x+mu).*y + 3*mu./r.^5.*(x-1+mu).*y;
    U_yy(:) = 1 - (1-mu)./d.^3 - mu./r.^3 + 3*(1-mu)./d.^5.*y.^2 + 3*mu./r.^5.*y.^2;
    U_zz(:) =- (1-mu)./d.^3 - mu./r.^3 + 3*(1-mu)./d.^5.*z.^2 + 3*mu./r.^5.*z.^2;
    U_xz(:) = 3*(1-mu)./d.^5.*(x+mu).*z + 3*mu./r.^5.*(x-1+mu).*z;
    U_yz(:) = 3*(1-mu)./d.^5.*z.*y + 3*mu./r.^5.*z.*y;


    
end