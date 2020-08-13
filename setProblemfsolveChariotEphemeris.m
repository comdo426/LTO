function [Problem, State] = setProblemfsolveChariotEphemeris(System, State, Spacecraft, Option, ...
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
%   Date: 01-Feb-2020; Last revision: 01-Mar-2020

nPhase = length(State);

%% set initial guess
Problem.x0 = [];
Body = Option.Body;
JDFix = Option.JDFix;

for iPhase = 1:nPhase
	
	lstar = System{iPhase}.parameter.lstar;
	tstar = System{iPhase}.parameter.tstar;	
	
	nSegment = State{iPhase}.nSegment;
	stateMat = State{iPhase}.state;
	
	controlMat = State{iPhase}.control;
	Tmax = Spacecraft{iPhase,1}.thrustMaxND;
	lambda = nan(nSegment, 1);	
		
	for iSegment = 1:nSegment
		T = controlMat(iSegment, 1);
		if T <= Tmax && T >= 0
			lambda(iSegment, 1) = asin(sqrt(T/Tmax));
		else
			if T < 0
				lambda(iSegment, 1) = 0;
			else
				lambda(iSegment, 1) = pi/2; % when lambda becomes imag., make it "1"
			end
		end % T <= Tmax if loop
		if imag(lambda(iSegment,1)) ~=0
			save('TEST3')
			error('lambda wrong')
		end
	end % iSegment for loop
	
	stateArray = reshape(stateMat', [7*4*nSegment, 1]);
	controlArray = reshape(controlMat', [4*nSegment, 1]);
	
	Problem.x0 = [Problem.x0; stateArray; controlArray; lambda];
	
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

% Add the aop and mf as the variable. 
Problem.x0 = [Problem.x0; Option.SpiralDown.aop; Option.SpiralDown.mf; Option.SpiralDown.time];

%% get Ephemeris state of the gravitational bodies

	
for iPhase = 1:nPhase

	nSegment = State{iPhase}.nSegment;
	tVarVec = State{iPhase}.timeVariable;
	tDefVec = State{iPhase}.timeDefect;
	
	lstar = System{iPhase}.parameter.lstar;
	tstar = System{iPhase}.parameter.tstar;	
	
	% Make sure that these are row vectors, not column vectors
	
	if iPhase == 1
		JDFixPhase = JDFix;
	else
		JDFixPhase = JDFixPhase + ...
			State{iPhase-1}.timeSegment(end)*System{iPhase-1}.parameter.tstar/60/60/24;
	end
	
	etVarVec = (JDFixPhase - cspice_j2000)*60*60*24 + ...
		tVarVec'*System{iPhase}.parameter.tstar;
	
	etDefVec = (JDFixPhase - cspice_j2000)*60*60*24 + ...
		tDefVec'*System{iPhase}.parameter.tstar;
	
	framePhase = Option.FrameSystem{iPhase}.frame;
	CBPhase = Option.FrameSystem{iPhase}.centralBody;
	
	for iBody = 1:length(Body.GM)
		stateVar = cspice_spkezr(Body.ID{iBody}, etVarVec, ...
			framePhase, 'NONE', CBPhase);	
		stateDef = cspice_spkezr(Body.ID{iBody}, etDefVec, ...
			framePhase, 'NONE', CBPhase);
		bodyStateMatVar{iBody} = reshape([stateVar(1:3,:)/lstar; ...
			stateVar(4:6,:)/lstar*tstar], [6, 4, nSegment]);
		bodyStateMatDef{iBody} = reshape([stateDef(1:3,:)/lstar; ...
			stateDef(4:6,:)/lstar*tstar], [6, 3, nSegment]);
	end
State{iPhase}.bodyStateMatVar = bodyStateMatVar;
State{iPhase}.bodyStateMatDef = bodyStateMatDef;
end

%% set objective function
Problem.objective = @(x) fsolveConstraintChariotEphemeris(x, System, State, Spacecraft, ...
	Option, Collocation);

	save('TESTChariotEphemeris');
%% set options
if ~isempty(Option.newton)
	Problem.options = Option.newton;
else
	Problem.options = Option.fsolve;
end

%% set solver
Problem.solver = 'fsolve';

end