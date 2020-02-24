function Problem = setProblemEphemerisCoast(stateMat, time, JDFix, iFix, ...
	Option)
%SETPROBLEMEPHEMERISCOAST - sets Problem structure for the fsolve/newtonRaphson
%in the ephemeris for the coast/ballistic motion
%
%  Syntax:
%     Problem = SETPROBLEMEPHEMERISCOAST(stateMat, time, JDFix, iFix, ...
%     System, Body, Option
%
%  Description:
%     Sets the problem structure for the fsolve/newtonRaphson in the ephemeris 
%		for the coast motion
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
%  See also: PERIODICORBIT2EPHEMERIS
%
%   Author: Beom Park
%   Date: 18-Feb-2020; Last revision: 18-Feb-2020

%% set initial guess
Problem.x0 = [];
nTime = length(time);
nSegment = nTime-1;
nState = 6*nSegment + 6;
Problem.x0 = nan(nState + nSegment + nTime,1);
stateVector = reshape(stateMat', [nState, 1]); % transpose to make it right
FrameSystem = Option.FrameSystem;
Body = Option.Body;

tIntegration = nan(nSegment,1);
for iSegment = 1:nSegment
	tIntegration(iSegment) = time(iSegment+1) - time(iSegment);
end

Problem.x0 = [stateVector; tIntegration; time];
N.time = nTime;
N.segment = nSegment;
N.state = nState;

%% set objective function
Problem.objective = @(x) fsolveConstraintEphemerisCoast(x, JDFix, iFix, ...
	FrameSystem, Body, N);

%% set options
if ~isempty(Option.newton)
	Problem.options = Option.newton;
else
	Problem.options = Option.fsolve;
end

%% set solver
Problem.solver = 'fsolve';

end