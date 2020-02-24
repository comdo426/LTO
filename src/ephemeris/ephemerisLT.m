function dYdt = ephemerisLT(t, Y, u, JD, Isp, g0, System, Body)
%EPHEMERISLT - computes the time derivative of the states in ephemeris with
%low-thrust
%
%  Syntax:
%     dYdt = EPHEMERISLT(t,y,JD, System, Body)
%
%  Description:
%     computes the time derivative of the states in ephemeris with low thrust
%
%  Inputs:
%		t - relative time w.r.t. the reference JD. in nondimensional unit defined
%		by System.tstar.
%		y - state in nondimentional unit defined by System.lstar, System.tstar,
%		System.centralBody. Spacecraft.
%			[7,1] = size(y)
%     u - control in nondimensional unit defined by System.lstar, System.tstar,
%     m0(initial spacecraft mass), thrust is in nondeimensional unit, alpha and
%     beta angles are difined with respect to the Equator/Ecliptic J2000 frame
%     depending on your frame of choice
%        [3,1] = size(u)
%		JD - Julian date of the reference
%		System - structure with characteristic parameters, mu, mu1, mu2, central
%		body, frame information
%		Body - structure with influential body names and corresponding
%		gravitational parameters%
%
%  Outputs:
%     dYdt - time derivative of the ephemeris state
%
%   Author: Beom Park
%   Date: 18-Feb-2020; Last revision: 18-Feb-2020

lstar = System.lstar;
tstar = System.tstar;
frame = System.frame;
CB = System.centralBody;
iCB = find(contains(Body.ID, CB), 1);
if isempty(iCB)
	error('cannot find central body in the Body structure')
end

nBody = length(Body.GM);

%% central body gravity
x = Y(1);
y = Y(2);
z = Y(3);
xdot = Y(4);
ydot = Y(5);
zdot = Y(6);
m = Y(7);

T = u(1);
alpha = u(2);
beta = u(3);

A1 = 0;
A2 = 0;
A3 = 0;

for iBody = 1:nBody
	
	muND = Body.GM(iBody)*tstar^2/lstar^3;
	ID = Body.ID{iBody};
	
	isCentral = (iBody == iCB);
	
	if isCentral
		
		r = sqrt(x^2 + y^2 + z^2);
		
		A1 = A1 - x*muND/r^3 + T./m*cos(alpha)*cos(beta);
		A2 = A2 - y*muND/r^3 + T./m*sin(alpha)*cos(beta);
		A3 = A3 - z*muND/r^3 + T./m*sin(beta);
		mdot = -T/(Isp*g0);
		
	else
		
		secPastJ2000 = (JD-cspice_j2000)*60*60*24 + t*tstar;
		bodyRelative2Central = cspice_spkezr(ID, secPastJ2000, frame, 'NONE', CB);
		bodyRelative2CentralND(1:3) = bodyRelative2Central(1:3)/lstar;
		xRel = bodyRelative2CentralND(1)-x;
		yRel = bodyRelative2CentralND(2)-y;
		zRel = bodyRelative2CentralND(3)-z;
		rRel = sqrt(xRel^2 + yRel^2 + zRel^2);
		
		xCentral2BodyND = bodyRelative2CentralND(1);
		yCentral2BodyND = bodyRelative2CentralND(2);
		zCentral2BodyND = bodyRelative2CentralND(3);
		rBodyRel = sqrt(xCentral2BodyND^2 + yCentral2BodyND^2 + ...
			zCentral2BodyND^2);
		
		A1 = A1 + xRel*muND/rRel^3 - xCentral2BodyND*muND/rBodyRel^3;
		A2 = A2 + yRel*muND/rRel^3 - yCentral2BodyND*muND/rBodyRel^3;
		A3 = A3 + zRel*muND/rRel^3 - zCentral2BodyND*muND/rBodyRel^3;
	end
	
end

dYdt = [xdot; ydot; zdot; A1; A2; A3; mdot];

end