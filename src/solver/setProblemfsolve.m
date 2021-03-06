function [Problem, State] = setProblemfsolve(System, State, Spacecraft, Option, ...
	Collocation)
%SETPROBLEMFSOLVE - sets Problem structure for the fsolve/newtonRaphson
%
%  Syntax:
%     [Problem, State] = SETPROBLEMFSOLVE(System, State, Spacecraft, Option, ...
%		Collocation)
%
%  Description:
%     Produces Problem structure for the fsolve/newtonRaphson. It's composed of
%     four main sections.
%			- initialization: generate x0(initial guess) for fsolve and
%			newtonRaphson. Note that the size differs if we have additional
%			constraints(e.g., altitude constraints)
%			- set objective function: sets the function to be solved for zero. In
%			this function it always calls fsolveConstraint.m to be its objective
%			function
%			- set options: pass on the fsolve/newtonRaphson option from the struct
%			Option
%			- set solver: this is set to be 'fsolve'. When using newtonRaphson,
%			this solver will not be needed
%
%  Outputs:
%     Problem - structure with following variables
%			x0 - column vector for the initial guess
%			options - structure for the fsolve/newtonRaphson option
%			solver - set to be 'fsolve'
%			objective - the function to solve for zero
%		State - after the slack variables are added
%
%  See also: LTOMAIN, SETPROBLEMFMINCON, FSOLVECONSTRAINTS
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 16-Feb-2020

nPhase = length(State);

%% set initial guess
Problem.x0 = [];

for iPhase = 1:nPhase
	[~, nSegment, ~, ~, ~, ~, ~] = getPhaseStateInfo(State{iPhase});
	Tmax = Spacecraft{iPhase,1}.thrustMaxND;
	lambda = nan(nSegment, 1);
	[~, ~, controlMat] = getStateControlMat(State{iPhase});
	for iSegment = 1:nSegment
		T = controlMat(iSegment, 1);
		if T <= Tmax
			lambda(iSegment, 1) = asin(sqrt(T/Tmax));
		else
			lambda(iSegmemt, 1) = pi/2; % when lambda becomes imag., make it "1"
		end % T <= Tmax if loop
	end % iSegment for loop
	Problem.x0 = [Problem.x0; State{iPhase}.state; ...
		State{iPhase}.control; lambda];
	State{iPhase}.slack = [];
	State{iPhase}.slack = lambda;
	
	% currently only supports altitude constraint
	if ~isempty(Option.AddCon{iPhase, 1})
		mnAddCon = size(Option.AddCon);
		nAddCon = mnAddCon(1, 2);
		for iAddCon = 1:nAddCon
			Con = Option.AddCon{iPhase, iAddCon};
			isSlack = strcmp(Con{2}, 'ineq'); % If slack, we need to augment x0
			if isSlack
				c = str2func(strcat('set',Con{3},'Constraint')); % define function
				sigma = c(System{iPhase}, State{iPhase}, Spacecraft{iPhase}, Con);
				Problem.x0 = [Problem.x0; sigma];
				State{iPhase}.slack = [State{iPhase}.slack; sigma];
			end
		end
	end % Additional constraint if loop
	
end % iPhase for loop

%% set objective function
Problem.objective = @(x) fsolveConstraint(x, System, State, Spacecraft, ...
	Option, Collocation);

%% set options
if ~isempty(Option.newton)
	Problem.options = Option.newton;
else
	Problem.options = Option.fsolve;
end

%% set solver
Problem.solver = 'fsolve';

end