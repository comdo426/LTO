function sigma = setAltitudeConstraint(SystemPhase, StatePhase, ~, Con)
%SETALTITUDECONSTRAINT - sets the slack variable for the altitude constraints
%
%  Syntax:
%     sigma = SETALTITUDECONSTRAINT(SystemPhase, StatePhase, ~, Con)
%
%  Description:
%     sets the slack variable for the altitude constraints from the bodies. Note
%     that the contraint is set to be x^2 + y^2 + z^2 >= h^2, to make the
%     partial computation easier
%
%  Inputs:
%		Con - cell that contains the altitude information
%			Con{4} - number of the primary(can be either 1/2 for CR3BP)
%			Con{5} - distance from the primary(NOT the altitude, the distance from
%			the center of the primary. You need to set the radius from the script
%			before you enter LTOMain
%
%  Outputs:
%     sigma - vector that contains the slack variable for the altitude
%     constraints
%
%  See also: SETPROBLEMFSOLVE, SETPROBLEMFMINCON, MULTIPHASESAMEDYNAMICS
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 16-Feb-2020

mu = SystemPhase.parameter.mu;
lstar = SystemPhase.parameter.lstar;

[~, ~, nNode, ~, ~, ~, ~] = getPhaseStateInfo(StatePhase);
stateMat = StatePhase.state;

primaryNo = Con{4};
minAlt = Con{5}/lstar;

isCR3BP = strcmp(SystemPhase.dynamics{1}, 'CR3BP');

sigma = nan(nNode, 1);

for iNode = 1:nNode
	if isCR3BP
		if primaryNo == 1
			alt = norm(stateMat(iNode, 1:3) - [-mu, 0, 0]);
		else
			if primaryNo == 2
				alt = norm(stateMat(iNode, 1:3) - [1-mu, 0, 0]);
			else
				error('primaryNo input Wrong');
			end
		end % primaryNo == 1 if loop
	else
		error('Dynamics other than CR3BP not supported');
	end % isCR3BP if loop
	
	if alt >= minAlt
		sigma(iNode, 1) = asin(minAlt^2/alt^2);
	else
		sigma(iNode, 1) = pi/2; % when sigma becomes imag., make it "1"
		fprintf('Altitude constraint violation at %d\n', iNode);
	end % alt >= minAlt if loop
end

end