function stateRev = propagateStateInitial(s, tSegment, stateInitial, ...
	includeLast, SystemPhase)
%PROPAGATESTATEINITIAL - propagates with a givn initial state(6D)
%
%  Syntax:
%     stateRev = PROPAGATESTATEINITIAL(s, tSegment, stateInitial, stateInitial,
%     includeLast, SystemPhase)
%
%  Description:
%     Produces a state matrix(each row is a 7D vector with position, velocity,
%     and the mass) at each node within the LGL 7th method collocation
%     algorithm.
%
%  Inputs:
%     s - number of segment over the given time
%			[1,1] = size(s); double = class(s)
%		tSegment - time for the segment
%			[s+1,1]  = size(tSegment); double = class(tSegment)
%		stateInitial - initial state for the propagation
%			[1,6] = size(stateInitial); double = class(stateInitial)
%		includeLast - boolean variable, indicating whether the stateRev matrix
%			should include the last propagation. For example, if includeLast == 1,
%			then stateRev(end, :) is the state at tSegment(end). Otherwise,
%			stateRev(end, :) is the state at tSegment(end-1).
%			[1,1] = size(includeLast); logical = class(includeLast)
%		SystemPhase - System at each phase. This is needed for the mu value of the
%			CR3BP propagation. 
%			[1,1] = size(SystemPhase); struct = class(SystemPhase)
%
%  Outputs:
%     stateRev - matrix with the propagated states 
%			[3*s+1, 7] OR [3*s, 7] = size(stateRev); double = class(stateRev)
%			Note that the size depends on the boolean value of includeLast
%
%
%  Other m-files required: CR3BP
%
%   Author: Beom Park
%   Date: 16-Feb-2020; Last revision: 16-Feb-2020


opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-17);
tau = [-1; -sqrt(495 +66*sqrt(15))/33; -sqrt(495 -66*sqrt(15))/33; 0; ...
	sqrt(495 -66*sqrt(15))/33; sqrt(495 +66*sqrt(15))/33; 1];

stateRev(1,:) = [stateInitial, 1];

mu = SystemPhase.parameter.mu;

for iSegment = 1:s
	tNode = zeros(7,1);
	dtSegment = tSegment(iSegment+1) - tSegment(iSegment);
	for itau = 1:7
		tNode(itau) = tSegment(iSegment) + dtSegment*(tau(itau)+1)/2;
	end % itau for loop
	[~, state] = ode113(@(t,y) CR3BP(t,y,mu,1,0), [tNode(1), tNode(3)], ...
		stateRev(3*(iSegment-1)+1, 1:6), opts);
	stateRev(3*(iSegment-1)+2, :) = [state(end, :), 1];
	[~, state] = ode113(@(t,y) CR3BP(t,y,mu,1,0), [tNode(3), tNode(5)], ...
		stateRev(3*(iSegment-1)+2, 1:6), opts);
	stateRev(3*(iSegment-1)+3, :) = [state(end, :), 1];
	if iSegment < s || includeLast
		[~, state] = ode113(@(t,y) CR3BP(t,y,mu,1,0), [tNode(1), tNode(3)], ...
			stateRev(3*(iSegment-1)+3, 1:6), opts);
		stateRev(3*(iSegment-1)+4, :) = [state(end, :), 1];
	end
end % iSegment for loop
end