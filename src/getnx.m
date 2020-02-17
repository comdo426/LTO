function [nx, nxOpt] = getnx(State)
%GETNX - gets the number of parameters for the fsolve/newtonRaphson/fmincon
%
%  Syntax:
%     [nx, nxOpt] = GETNX(State)
%
%  Description:
%     gets the number of parameters for the each iPhase. Note that the nx, nxOpt
%     will give you different numbers because State{iPhase}.slack is defined
%     separately from setProblem functions. 
%
%  Outputs:
%     nx - column vector of numbers of parameters for each phase
%		nxOpt - column vector of numbers of parameters for each phase(optimizer)
%
%  See also: SETPROBLEMFSOLVE, SETPROBLEMFMINCON
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 16-Feb-2020

nPhase = length(State);
nx = zeros(nPhase, 1);
nxOpt = zeros(nPhase, 1);

for iPhase = 1:nPhase
	[~, nSegment, ~, nState, nControl, ~, nSlack] = getPhaseStateInfo(State{iPhase});
	nx(iPhase) = nState + nControl + nSlack;
	nxOpt(iPhase) = nState + nControl + nSlack;
	% works for altitude constraints
end

end