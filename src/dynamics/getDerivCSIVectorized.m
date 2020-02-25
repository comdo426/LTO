function Ydot = getDerivCSIVectorized(mu, Y, U, isp, g0)
%GETDERIVCSIVECTORIZED - computes the F and dF matrix for fsolve/newtonRaphson
%
%  Syntax:
%     Ydot = GETDERIVCSIVECTORIZED(mu, Y, U, alpha, beta, isp, g0)
%
%  Description:
%     Computes derivative of CSI, in the vectorized form. The inputs should be
%     all vectorized, and are 3D arrays. 
%
%  Inputs:
%		Y - state.
%        [7, 3~4, nSegment] = size(Y)
%        size(Y,1) - 7, since we are looking at the position, velocity, mass
%        size(Y,2) - 3 if we're computing defect nodes, 4 if we're computing
%        variable nodes
%        size(Y,3) - nSegment, corresponds to the different pages
%     U - control
%        [3, 1, nSegment] = size(U)
%     isp - nondimensional specific impulse
%     g0 - nondimensional gravitational acceleration at Earth ground
%       
%  Outputs:
%     Ydot - state derivative w.r.t. time
%        [7, 3~4, nSegment] = size(Ydot)
%
%  See also: FMINCONCONSTRAINT
%
%   Author: Beom Park
%   Date: 24-Feb-2020; Last revision: 24-Feb-2020


x = Y(1, :, :);
y = Y(2, :, :);
z = Y(3, :, :);
xdot = Y(4, :, :);
ydot = Y(5, :, :);
zdot = Y(6, :, :);
m = Y(7, :, :);

T = U(1, :, :);
alpha = U(2, :, :);
beta = U(3, :, :);

d = sqrt((x+mu).^2 + y.^2 + z.^2);
r = sqrt((x -1 + mu).^2 +y.^2 + z.^2);

A1 = 2.*ydot + x - (1-mu)./d.^3.*(x+mu) - mu./r.^3.*(x-1+mu) + T./m.*cos(alpha).*cos(beta);
A2 = -2.*xdot + y - (1-mu)./d.^3.*y - mu./r.^3.*y + T./m.*sin(alpha).*cos(beta);
A3 = - (1-mu)./d.^3.*z - mu./r.^3.*z + T./m.*sin(beta);

mdot = -T/isp/g0;

Ydot = [xdot; ydot; zdot; A1; A2; A3; mdot];

end