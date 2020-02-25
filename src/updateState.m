function State = updateState(State, x)
%UPDATESTATE - update the converged solution to the State cell
%
%  Syntax:
%     State = UPDATESTATE(State, x)
%
%  Description:
%     replace the vectors within State cell to the converged solution
%
%  Inputs:
%		State - old State cell to be updated
%		x - converged solution from feasible/optimize solvers
%
%  Outputs:
%     State - new State cell after the update
%
%  See also: GETPHASESTATEINFO, GETNX
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 16-Feb-2020

nPhase = length(State);

% reset slack vectors to blank
for iPhase = 1:nPhase
	State{iPhase}.slack = [];
end

nx = getnx(State);

for iPhase = 1:nPhase
	[~, nSegment, ~, nState, nControl, ~, ~] = getPhaseStateInfo(State{iPhase});
	
	if iPhase == 1
		nb = 0;
	else
		nb = sum(nx(1:iPhase-1));
	end
	State{iPhase}.state = reshape(x(1+nb:nState+nb), [7, 4*nSegment])';
	State{iPhase}.control = reshape(x(1+nState+nb:nState+nControl+nb), [3, nSegment])';
end % iPhase for loop

end