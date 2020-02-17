function [Problem, State] = ...
	setProblemfmincon(System, State, Spacecraft, Option, Collocation)
%SETPROBLEMFMINCON - sets Problem structure for fmincon
%
%  Syntax:
%     [Problem, State] = ...
%		SETPROBLEMFMINCON(System, State, Spacecraft, Option, Collocation)
%
%  Description:
%     Produces Problem structure for the fmincon. It's composed of
%     following sections.
%			- initialization: generate x0(initial guess) for fmincon. Note that the
%			size differs depending on the number of constraints
%			constraints(e.g., altitude constraints)
%			- set objective function: sets the function to be minimized. In
%			this function it always calls fminconObjective.m
%			- set linear constraint(optional)
%			- set boundary constraint(optional)
%			- set nonlinear constraint function. fminconConstraint function called.
%			- set option: pass on the user input Option
%			- set solver: this is set to be 'fmincon'
%
%  Outputs:
%     Problem - structure with properties for the Problem
%		State - after the slack variables are added
%
%  See also: SETPROBLEMFSOLVE, FMINCONOBJECTIVE, FMINCONCONSTRAINT
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 16-Feb-2020

nPhase = length(State);

%% initialize

Problem.x0 = [];

for iPhase = 1:nPhase
	
	Problem.x0 = [Problem.x0; State{iPhase}.state; State{iPhase}.control];
	State{iPhase}.slack = [];
	
	if ~isempty(Option.AddCon{iPhase, 1})
		mnAddCon = size(Option.AddCon);
		nAddCon = mnAddCon(1, 2);
		
		for iAddCon = 1:nAddCon
			Con = Option.AddCon{iPhase, iAddCon};
			isNlnr = strcmp(Con{1}, 'nlnr');
			isSlack = strcmp(Con{2}, 'ineq'); % If slack, we need to augment the x0
			if isNlnr && isSlack
				c = str2func(strcat('set',Con{3},'Constraint'));
				sigma = c(System{iPhase}, State{iPhase}, Spacecraft{iPhase}, Con);
				Problem.x0 = [Problem.x0; sigma];
				State{iPhase}.slack = [State{iPhase}.slack; sigma];
			end
		end
	end
	
end

[~, nxOpt] = getnx(State);

iFinalMass = 0;

for iPhase = 1:nPhase
	[~, ~, ~, nState, ~, ~, ~] = getPhaseStateInfo(State{iPhase});
	if iPhase < nPhase
		iFinalMass = iFinalMass + nxOpt(iPhase);
	else
		iFinalMass = iFinalMass + nState;
	end
	
end


%% set objective function
Problem.objective = @(x) fminconObjective(x, iFinalMass);

%% set linear equality constraints

% These were commented, since I am using m = 1 as the continuity constraint
% inside the fminconConstraint. If this were to change, make sure you uncomment
% these lines

% xLength = length(Problem.x0);
% Aeq = zeros(1, xLength);
% Aeq(1, 7) = 1; % indicates the initial mass
% Problem.Aeq = Aeq;
% Problem.beq = 1; % initial mass should be 1

%% set boundary for parameters

lb = [];
ub = [];

for iPhase = 1:nPhase
	[~, nSegment, ~, nState, ~, ~, ~] = getPhaseStateInfo(State{iPhase});
	lbPhase = -inf*ones(nxOpt(iPhase), 1);
	ubPhase = inf*ones(nxOpt(iPhase), 1);
	for iSegment = 1:nSegment
		lbPhase(nState + 3*(iSegment-1) + 1) = -eps;
		ubPhase(nState + 3*(iSegment-1) + 1) = ...
			Spacecraft{iPhase}.thrustMaxND+eps;
	end
	lb = [lb; lbPhase];
	ub = [ub; ubPhase];
	clear lbPhase ubPhase
end

Problem.lb = lb;
Problem.ub = ub;

%% set nonlinear constraint

Problem.nonlcon = @(x) fminconConstraint(x, System, State, Spacecraft, ...
	Option, Collocation); 

%% set options
Problem.options = Option.fmincon;
%% set solver
Problem.solver = 'fmincon';
end