function [stateInertConvergedND, timeConverged] = ...
	periodicOrbit2Ephemeris(SystemPhase, state, ...
	period, Option, EphemOptionPhase)
%PERIODICORBIT2EPHEMERIS - converge periodic orbits in ephemeris
%
%  Syntax:
%     [stateInertConvergedND, timeConverged] = ...
%	PERIODICORBIT2EPEHEMERIS(SystemPhase, state, ...
%	period, Option, EphemOptionPhase)
%
%  Description:
%     Stack periodic orbits by nRev / s segments and converge them in ephemeris
%
%  Outputs:
%		stateInertConvergedND - state converged in the ephemeris, in the frame
%		set in Option.FrameSystem. They are non-dimensionalized in the
%		charateristic quantities defined in Option.FrameSystem
%		timeConverged - time converged in nondimensional unit w.r.t.
%		Option.FrameSystem
%
%  See also: PERIODICORBIT2EPHEMERIS
%
%   Author: Beom Park
%   Date: 17-Feb-2020; Last revision: 17-Feb-2020

mu = SystemPhase.parameter.mu;
lstar = SystemPhase.parameter.lstar;
tstar = SystemPhase.parameter.tstar;

nRev = EphemOptionPhase.nRev;
s = EphemOptionPhase.s;
JDFix = EphemOptionPhase.JD;
tOffset = EphemOptionPhase.tOffset;
FrameSystem = Option.FrameSystem;
opts = Option.integrate;
iFix = floor(nRev/2)*s + 1; % index of the state to be fixed along the stack
tOneRev = linspace(0, period, s+1)';

% Stack the ballistic CR3BP orbits

if tOffset ~= 0
	[~, Y] = ode113(@(t,y) CR3BP(t,y,mu,1,0), ...
		[0, tOffset], state, opts);
	stateOneRev(1,:) = Y(end,:);
else
	stateOneRev(1,:) = state';
end

for i = 1:s
	[~, Y] = ode113(@(t,y) CR3BP(t,y,mu,1,0), ...
		[tOneRev(i), tOneRev(i+1)], stateOneRev(i,:), opts);
	stateOneRev(i+1,:) = Y(end,:);
end

stateRotTTL = [];
timeTTL = [];

for iRev = 1:nRev
	stateRotTTL = [stateRotTTL; stateOneRev(1:end-1,:)]; %#ok<AGROW>
	timeTTL = [timeTTL; tOneRev(1:end-1)+period*(iRev-1)]; %#ok<AGROW>
end

% Add the final state to complete the final revolution
stateRotTTL = [stateRotTTL; stateOneRev(1, :)];
timeTTL = [timeTTL; period*nRev];

% move the time w.r.t. the time to be fixed
timeTTL = timeTTL - timeTTL(iFix);
timeTTLSecond = timeTTL*tstar;

% state in inertial frame, dimensional units
stateInert = rot2inert(stateRotTTL, timeTTLSecond, JDFix, FrameSystem);
% uncomment this to check visually
% figure(3)
% plot3(stateInert(:,1), stateInert(:,2), stateInert(:,3), 'k')

% non-dimensionalization
stateInertND = [stateInert(:,1:3)/lstar, stateInert(:,4:6)/lstar*tstar];
Problem = setProblemEphemerisCoast(stateInertND, timeTTL, JDFix, ...
	iFix, Option);

% Choose the solver
isNewton = ~isempty(Option.newton);
if isNewton
	xConverged = newtonRaphson(Problem);
else
	xConverged = fsolve(Problem);
end

% parse the converged vector back into matrix form
nTime = length(timeTTL);
nSegment = nTime-1;
nState = 6*nSegment+6;
stateConvergedVec = xConverged(1:nState);
stateInertConvergedND = reshape(stateConvergedVec, [6, nTime])'; ...
	%Transpose to match the dimensions
timeConverged = xConverged(nState+nSegment+1:nState+nSegment+nTime);

% Uncomment the following lines to check the states in the inertial frame 

% stateInertConverged = [stateInertConvergedND(:,1:3)*lstar, ...
% 	stateInertConvergedND(:,4:6)*lstar/tstar];
% figure(2)
% plot3(stateInertConverged(:,1), stateInertConverged(:,2), stateInertConverged(:,3))

end

