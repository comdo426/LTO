function Ydot = getDerivCSI(mu, Y, T, alpha, beta, isp, g0)

x = Y(1, :);
y = Y(2, :);
z = Y(3, :);
xdot = Y(4, :);
ydot = Y(5, :);
zdot = Y(6, :);
m = Y(7, :);

d = sqrt((x+mu).^2 + y.^2 + z.^2);
r = sqrt((x -1 + mu).^2 +y.^2 + z.^2);

A1 = 2.*ydot + x - (1-mu)./d.^3.*(x+mu) - mu./r.^3.*(x-1+mu) + T./m.*cos(alpha).*cos(beta);
A2 = -2.*xdot + y - (1-mu)./d.^3.*y - mu./r.^3.*y + T./m.*sin(alpha).*cos(beta);
A3 = - (1-mu)./d.^3.*z - mu./r.^3.*z + T./m.*sin(beta);

mdot = -T/isp/g0;

Ydot = [xdot; ydot; zdot; A1; A2; A3; mdot];

end