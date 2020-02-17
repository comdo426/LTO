function xShrink = deleteSlackVariable(x, State, ~)
%DELETESLACKVARIABLE - NewtonRaphson method to solve for zero
%
%  Syntax:
%     xShrink = DELETESLACKVARIABLE(x, State, Option, isFeasible, isOptimize)
%
%  Description:
%     Uses NewtonRaphson method to solve for zero of the objective function.
%     Note that it used to be different for feasible/optimize solvers, but I
%     changed them to be simple. 
%
%  Inputs:
%		x - converged vector to be shrinked
%		isFeasible/isOptimize - boolean
%
%  Outputs:
%     xShrink - vector after its slack variables are all deleted.
%
%  See also: SETPROBLEMFSOLVE, SETPROBLEMFMINCON
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 16-Feb-2020

nPhase = length(State);

nx = getnx(State);
nb = nan(nPhase, 1);


% compute nb: number of parameters before the phase
for iPhase = 1:nPhase
	if iPhase == 1
		nb(1) = 0;
	else
		nb(iPhase) = sum(nx(1:iPhase-1));
	end
end

for iPhase = 1:nPhase
	[~, ~, ~, nState, nControl, ~, nSlack] = ...
		getPhaseStateInfo(State{iPhase});
	x(nState+nControl+1+nb(iPhase) : nState+nControl+nSlack+nb(iPhase)) = [];
	nb = nb-nSlack;
end

xShrink = x;
end