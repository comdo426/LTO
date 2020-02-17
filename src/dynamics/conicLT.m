function dYdt = conicLT(~, Y, u1, ~, Isp, g0, ~)
%CONICLT - CSI, Inertial propagator for 2BP
%
%  Syntax:
%     dYdt = CONICLT(~, Y, u1, ~, Isp, g0, ~)
%
%  Description:
%     CSI, Inertial propagator for 2BP. Note that some of the input variables
%     are missing, as it has the same structure with CR3BPLT but has less
%     inputs.
%
%  Inputs:
%		Y - state(n.d.) consisted of 3 position, 3 velocity and 1 mass
%		u1 - control variables(n.d.) consisted of thrust, and alpha/beta(angles).
%		Note that angles are here defined in the J2000 equator coordinates
%		Isp - specific impulse(n.d.)
%		g0 - Earth gravitational acceleration at ground level(n.d.)
%
%  Outputs:
%     dYdt - time derivative of states
%
%	See also: CR3BPLT
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
    alpha = u1(2);
    beta = u1(3);
	 
	 r = sqrt(x.^2 + y.^2 + z.^2);
	 
	 A1 = -x/r^3 + T./m.*cos(alpha).*cos(beta);
	 A2 = -y/r^3 + T./m.*sin(alpha).*cos(beta);
	 A3 = -z/r^3 + T./m.*sin(beta);
	 
	 mdot = -T/Isp/g0;
	 
	 dYdt = [xdot; ydot; zdot; A1; A2; A3; mdot];


end