function [x0, funcs, ipoptopt, State] = setProblemipoptEphemeris(System, State, ...
	Spacecraft, Option, Collocation)
%SETPROBLEMIPOPTEPHEMERIS- sets Problem structure for ipopt
%in the ephemeris with the low-thrust capabilities using collocation
%
%  Syntax:
%     Problem = SETPROBLEMIPOPTEPHEMERIS(stateMat, time, JDFix, iFix, ...
%     System, Body, Option
%
%  Description:
%     Sets the problem structure for the fsolve/newtonRaphson in the ephemeris
%		with the low-thrust capabilities
%
%	Inputs:
%		state - matrix that contains the inertial state with respect to the
%		central body defined as System.centralBody. Since this is a coasting
%		problem, mass is NOT included. All the units should be non-dimensional
%		with respect to System.lstar, System.tstar
%			[nSegment+1,6] = size(state)
%		time - colmun vector that contains the relative time from JDFix. The unit
%		is non-dimensional with respect to System.tstar
%			[nSegmetn,1] = size(time)
%		JDFix - Julian date to be fixed(days)
%		iFix - number of the state to be fixed(count)
%		System - structure that contains following information. Note that the
%		following only corresponds to CR3BP notations
%			P1 - string of the first primary
%			P2 - string of the second primary
%			mu - mass ratio
%			lstar - characteristic length (km)
%			tstar - characteristic length (km)
%			mu1 - gravitational parameter for the first primary
%			mu2 - gravitational parameter for the second primary
%			centralBody - centralBody
%			frame - string of the central body of motion
%		Body - structure that participates in the gravity
%			GM - gravitational parameter (km^2/sec^3)
%			ID - name of the bodies (cell)
%
%  Outputs:STATE
%     Problem - structure with following variables
%			x0 - column vector for the initial guess
%			options - structure for the fsolve/newtonRaphson option
%			solver - set to be 'fsolve'
%			objective - the function to solve for zero
%
%  See also: PERIODICORBIT2EPHEMERIS, SETPROBLEMEPHEMERISCOAST
%
%   Author: Beom Park
%   Date: 02-Mar-2020; Last revision: 02-Mar-2020

%% set initial guess

Problem.x0 = [];
nPhase = length(System);
JDFix = Option.JDFix;
stateFix = Option.stateFix;
FrameSystem = Option.FrameSystem;
lstar = FrameSystem.lstar;
tstar = FrameSystem.tstar;
Body = Option.Body;
C = Option.C;

x0 = [];
ipoptopt = Option.ipoptopt;

for iPhase = 1:nPhase
	nSegment = State{iPhase}.nSegment;
	stateMat = State{iPhase}.state;
	controlMat = State{iPhase}.control;
	stateArray = reshape(stateMat', [7*4*nSegment, 1]);
	controlArray = reshape(controlMat', [4*nSegment, 1]);
	x0 = [x0; stateArray; controlArray];
	State{iPhase}.slack = [];
		
	
	% Altitude constraint currently NOT supported
end

% 
% save('TEST2')
% error('TEST')

%% Get the index for the final mass
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



%% set boundary for parameters

lb = [];
ub = [];

posBoundaryD = Option.posBoundary; % km;
velBoundaryD = Option.velBoundary; % km/s;

for iPhase = 1:nPhase
	[~, nSegment, nNode, nState, ~, ~, ~] = getPhaseStateInfo(State{iPhase});
	stateMat = State{iPhase}.state;
	Tmax = Spacecraft{iPhase,1}.thrustMaxND;
	
	lbPhase = -inf*ones(nxOpt(iPhase), 1);
	ubPhase = inf*ones(nxOpt(iPhase), 1);
	
	% Non-dimensionalize the relative bounds set
	posBoundaryND = posBoundaryD/lstar;
	velBoundaryND = velBoundaryD/lstar*tstar;
		
	for iNode = 1:nNode
		
		% position, velocity, mass bounds
		lbPhase(7*(iNode-1)+1:7*(iNode-1)+3) = stateMat(iNode, 1:3) - posBoundaryND*ones(1,3);
		lbPhase(7*(iNode-1)+4:7*(iNode-1)+6) = stateMat(iNode, 4:6) - velBoundaryND*ones(1,3);
		lbPhase(7*(iNode-1)+7) = 0;
		
		ubPhase(7*(iNode-1)+1:7*(iNode-1)+3) = stateMat(iNode, 1:3) + posBoundaryND*ones(1,3);
		ubPhase(7*(iNode-1)+4:7*(iNode-1)+6) = stateMat(iNode, 4:6) + velBoundaryND*ones(1,3);
		ubPhase(7*(iNode-1)+7) = 1;
		
	end
	for iSegment = 1:nSegment
		
		% Control variables bounds
		lbPhase(nState + 4*(iSegment-1) + 1) = 0;
		ubPhase(nState + 4*(iSegment-1) + 1) = Tmax;
		
		lbPhase(nState + 4*(iSegment-1) + 2: nState + 4*(iSegment-1) + 4) = -ones(1,3);
		ubPhase(nState + 4*(iSegment-1) + 2: nState + 4*(iSegment-1) + 4) = ones(1,3);
		
	end
	
	lb = [lb; lbPhase];
	ub = [ub; ubPhase];
	clear lbPhase ubPhase
end

mContinuityTTL = 0;
isIniCon = ~isempty(State{1}.initialConstraint);
if isIniCon
	stateIniCon = State{1}.initialConstraint;
	nIniCon = length(stateIniCon);
	mContinuityTTL = mContinuityTTL + nIniCon;
end
isFinCon = ~isempty(State{end}.finalConstraint);
if isFinCon
	stateFinCon = State{end}.finalConstraint;
	nFinCon = length(stateFinCon);
	mContinuityTTL = mContinuityTTL + nFinCon;
end
mContinuityTTL = mContinuityTTL + 7*(nPhase-1);

[~, mcOpt] = getmc(System, State, Spacecraft, Option);
m = sum(mcOpt) + mContinuityTTL; % row number of constraint

% Define the bounds for the variables and the constraints
ipoptopt.lb = lb;
ipoptopt.ub = ub;
ipoptopt.cl = zeros(m,1);
ipoptopt.cu = zeros(m,1);

%% Define the auxdata structure for ipopt
ipoptopt.auxdata.System = System;
ipoptopt.auxdata.State = State;
ipoptopt.auxdata.Spacecraft = Spacecraft;
ipoptopt.auxdata.Option = Option;
ipoptopt.auxdata.Collocation = Collocation;
ipoptopt.auxdata.iFinalMass = iFinalMass;


%% set functions

funcs.objective         = @ipoptObjective;
funcs.constraints       = @ipoptConstraintEphemeris;
funcs.gradient          = @ipoptGradient;
funcs.jacobian          = @ipoptJacobianEphemeris;
funcs.jacobianstructure = @ipoptJacobianStructureEphemeris;

end