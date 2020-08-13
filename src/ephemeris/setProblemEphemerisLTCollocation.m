function [Problem, State] = setProblemEphemerisLTCollocation(System, State, ...
	Spacecraft, Option, Collocation)
%SETPROBLEMEPHEMERISLTCOLLOCATION - sets Problem structure for the fsolve/newtonRaphson
%in the ephemeris with the low-thrust capabilities using collocation
%
%  Syntax:
%     Problem = SETPROBLEMEPHEMERISCOAST(stateMat, time, JDFix, iFix, ...
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

for iPhase = 1:nPhase
	nSegment = State{iPhase}.nSegment;
	stateMat = State{iPhase}.state;
	tVar = State{iPhase}.timeVariable;
	
	posvelRot = stateMat(:, 1:6);
	mass = stateMat(:, 7);
		
	posvelInertD = rot2inert(posvelRot, tVar*tstar, ...
		JDFix.initial, FrameSystem);
	
	figure(150)
	hold on
	axis equal
	plot3(posvelInertD(:,1), posvelInertD(:,2), posvelInertD(:,3), 'r')
	
	posvelInertND = [posvelInertD(:, 1:3)/lstar, posvelInertD(:, 4:6)/lstar*tstar];
	stateInert = [posvelInertND, mass];
		
	controlMat = State{iPhase}.control;
	for i = 1:nSegment
		controlMat(i, 2:3) = (C*controlMat(i,2:3)')';
	end
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
	
	save('TEST')
	
	stateArray = reshape(stateInert', [7*4*nSegment, 1]);
	controlArray = reshape(controlMat', [4*nSegment, 1]);
	
	Problem.x0 = [Problem.x0; stateArray; controlArray; lambda];
	
	State{iPhase}.slack = lambda;
	
	if iPhase == 1
		State{iPhase}.initialConstraint = stateFix.initial;
	end
	if iPhase == nPhase
		State{nPhase}.finalConstraint = stateFix.final;
	end
	
	% Altitude constraint currently NOT supported
end

% 
% save('TEST2')
% error('TEST')


%% get Ephemeris state of the gravitational bodies

for iPhase = 1:nPhase

	nSegment = State{iPhase}.nSegment;
	tVarVec = State{iPhase}.timeVariable;
	tDefVec = State{iPhase}.timeDefect;
	
	% Make sure that these are row vectors, not column vectors
	etVarVec = (JDFix.initial - cspice_j2000)*60*60*24 + ...
		tVarVec'*System{iPhase}.parameter.tstar;
	
	etDefVec = (JDFix.initial - cspice_j2000)*60*60*24 + ...
		tDefVec'*System{iPhase}.parameter.tstar;
	
	for iBody = 1:length(Body.GM)
		stateVar = cspice_spkezr(Body.ID{iBody}, etVarVec, ...
			FrameSystem.frame, 'NONE', FrameSystem.centralBody);	
		stateDef = cspice_spkezr(Body.ID{iBody}, etDefVec, ...
			FrameSystem.frame, 'NONE', FrameSystem.centralBody);
		bodyStateMatVar{iBody} = reshape([stateVar(1:3,:)/lstar; ...
			stateVar(4:6,:)/lstar*tstar], [6, 4, nSegment]);
		bodyStateMatDef{iBody} = reshape([stateDef(1:3,:)/lstar; ...
			stateDef(4:6,:)/lstar*tstar], [6, 3, nSegment]);
	end
State{iPhase}.bodyStateMatVar = bodyStateMatVar;
State{iPhase}.bodyStateMatDef = bodyStateMatDef;
end


%% set objective function
Problem.objective =  @(x) fsolveConstraintEphemerisLTCollocation(x, System, ...
	State, Spacecraft, Option, Collocation);

%% set options
if ~isempty(Option.newton)
	Problem.options = Option.newton;
else
	Problem.options = Option.fsolve;
end

%% set solver
Problem.solver = 'fsolve';

end