function dYdt = CR3BPLT(t, Y, u1, mu, Isp, g0, t0)
%CR3BPLT - CSI, Inertial propagator for CR3BP + low-thrust
%
%  Syntax:
%     dYdt = CR3BPLT(t, Y, u1, mu, Isp, g0, t0)
%
%  Description:
%     CSI, Inertial propagator for CR3BP + low-thrust. 
%
%  Inputs:
%		t - time(n.d.)
%		Y - state(n.d.) consisted of 3 position, 3 velocity and 1 mass
%		u1 - control variables(n.d.) consisted of thrust, and alpha/beta(angles).
%		Note that angles are here defined in the J2000 equator coordinates
%		mu - mass ratio(n.d.)
%		Isp - specific impulse(n.d.)
%		g0 - Earth gravitational acceleration at ground level(n.d.)
%		t0 - time/angle offset(n.d.) to consider the phasing.
%
%  Outputs:
%     dYdt - time derivative of states
%
%	See also: CONICLT
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 17-Feb-2020

    x = Y(1);
    y = Y(2);
    z = Y(3);
    xdot = Y(4);
    ydot = Y(5);
    zdot = Y(6);
    m = Y(7);
    T = u1(1);
    alphaInert = u1(2);
    alpha = alphaInert - (t+t0); % subtract time, since angles are in inertial frame
    beta = u1(3);
    
    d = sqrt((x+mu).^2 + y.^2 + z.^2);
    r = sqrt((x -1 + mu).^2 +y.^2 + z.^2);
	 
    A1 = 2.*ydot + x - (1-mu)./d.^3.*(x+mu) - mu./r.^3.*(x-1+mu) + T./m.*cos(alpha).*cos(beta);
    A2 = -2.*xdot + y - (1-mu)./d.^3.*y - mu./r.^3.*y + T./m.*sin(alpha).*cos(beta);
    A3 = - (1-mu)./d.^3.*z - mu./r.^3.*z + T./m.*sin(beta);
        
    Vel = [xdot; ydot; zdot];
    Acc = [A1; A2; A3];
    mdot = -T/Isp/g0;
    
    dYdt = [Vel; Acc; mdot];


end

