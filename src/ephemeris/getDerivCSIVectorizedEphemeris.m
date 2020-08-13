function Ydot = getDerivCSIVectorizedEphemeris(muND, Y, U,  isp, g0, ...
	bodyState, iCB)
%GETDERIVCSIVECTORIZEDEPHEMERIS - get derivative of collocation at
%variable/defect nodes. For ephemeris
%
%  Syntax:
%     Ydot = GETDERIVCSIVECTORIZED(mu, Y, U, alpha, beta, isp, g0)
%
%  Description:
%     Computes derivative of CSI, in the vectorized form. The inputs should be
%     all vectorized, and are 3D arrays. 
%
%  Inputs:
%		muND - array of mu(nondimensional unit) participating in the gravity
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
%		bodyStatePhase - cell that contains the state of gravitational bodies
%       
%  Outputs:
%     Ydot - state derivative w.r.t. time
%        [7, 3~4, nSegment] = size(Ydot)
%
%  See also: FMINCONCONSTRAINT
%
%   Author: Beom Park
%   Date: 24-Feb-2020; Last revision: 24-Feb-2020

sizeY = size(Y);

x = Y(1, :, :);
y = Y(2, :, :);
z = Y(3, :, :);
xdot = Y(4, :, :);
ydot = Y(5, :, :);
zdot = Y(6, :, :);
m = Y(7, :, :);

T = U(1, :, :);
ux = U(2, :, :);
uy = U(3, :, :);
uz = U(4, :, :);

A1 = zeros(1, sizeY(2), sizeY(3));
A2 = zeros(1, sizeY(2), sizeY(3));
A3 = zeros(1, sizeY(2), sizeY(3));


for iBody = 1:length(bodyState)
	
	isCentral = (iBody == iCB);
	
	if isCentral
				
		r = sqrt(x.^2 + y.^2 + z.^2);
		A1 = A1 - x.*muND(iBody)./r.^3 + T./m.*ux;
		A2 = A2 - y.*muND(iBody)./r.^3 + T./m.*uy;
		A3 = A3 - z.*muND(iBody)./r.^3 + T./m.*uz;
		mdot = -T/isp/g0;
		
	else
		
		xRel2CB = bodyState{iBody}(1, :, :);
		yRel2CB = bodyState{iBody}(2, :, :);
		zRel2CB = bodyState{iBody}(3, :, :);
		rRel2CB = sqrt(xRel2CB.^2 + yRel2CB.^2 + zRel2CB.^2);
		
		xRel = bodyState{iBody}(1, :, :) - x;
		yRel = bodyState{iBody}(2, :, :) - y;
		zRel = bodyState{iBody}(3, :, :) - z;
		rRel = sqrt(xRel.^2 + yRel.^2 + zRel.^2);
		
		A1 = A1 + xRel.*muND(iBody)./rRel.^3 - xRel2CB.*muND(iBody)./rRel2CB.^3;
		A2 = A2 + yRel.*muND(iBody)./rRel.^3 - yRel2CB.*muND(iBody)./rRel2CB.^3;
		A3 = A3 + zRel.*muND(iBody)./rRel.^3 - zRel2CB.*muND(iBody)./rRel2CB.^3;
		

	end
	
end
% 
% save('TEST2')
% error('TEST')



Ydot = [xdot; ydot; zdot; A1; A2; A3; mdot];

end