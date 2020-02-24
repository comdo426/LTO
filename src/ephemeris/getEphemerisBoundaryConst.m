function [stateRotPlot,	JDFix, stateFix] = ...
	getEphemerisBoundaryConst(System, State, Option, EphemOption)
%GETEPHEMERISBOUNDARYCONST - get ephemeris boundary constraints
%
%  Syntax:
%     ephemerisBoundary = GETEPHEMERISBOUNDARYCONST(System, State, EphemOption)
%
%  Description:
%     Stack periodic orbits by nRev / s segments and converge them in ephemeris
%     to get the boundary constraints in ephemeris
%
%  Outputs:
%     stateRotPlot - structure that contains states in the rotating frame that
%     is propagated in ephemeris then rotated to the rotating frame. 
%		JDFix - structure that contains the fixed Julian Date
%		stateFix - structure that constains the fixed state
%
%  See also: PERIODICORBIT2EPHEMERIS
%
%   Author: Beom Park
%   Date: 17-Feb-2020; Last revision: 23-Feb-2020

%% initial periodic orbit

initialState = State{1}.initialConstraint(1:6);
initialPeriod = State{1}.period;
[stateInitialND, timeInitial] = periodicOrbit2Ephemeris(System{1}, ...
	initialState, initialPeriod, Option, EphemOption{1});

nSegment = length(timeInitial)-1;
stateInitialRotPlot = [];
JDFix.initial = EphemOption{1}.JD;
initialnRev = EphemOption{1}.nRev;
initials = EphemOption{1}.s;
% index of the state to be fixed along the stack
initialiFix = floor(initialnRev/2)*initials + 1; 

opts = Option.integrate;
FrameSystem = Option.FrameSystem;
lstar = FrameSystem.lstar;
tstar = FrameSystem.tstar;
Body = Option.Body;

% propagate the initial state to stack the ballistic CR3BP orbits
for i = 1:nSegment
	ts = [timeInitial(i), timeInitial(i+1)];
	Y0 = stateInitialND(i, :);
	
	[T, Y] = ode113(@(t,y) ephemeris(t,y,JDFix.initial, FrameSystem, Body), ts, Y0, opts);
	YRot = inert2rot([Y(:,1:3)*lstar, Y(:,4:6)*lstar/tstar], T*tstar, JDFix.initial, FrameSystem);
	stateInitialRotPlot = [stateInitialRotPlot; YRot]; %#ok<AGROW>
	if i == initialiFix
		stateFix.initial = [Y0, 1]';
		figure(Option.plotNo(1))
		hold on
		axis equal
		grid on
		plot3(YRot(1,1), YRot(1,2), YRot(1,3), '^r', 'linewidth', 2.0)
	end
	
end

figure(Option.plotNo(1))
plot3(stateInitialRotPlot(:,1), stateInitialRotPlot(:,2), stateInitialRotPlot(:,3), 'k')
stateRotPlot.initial = stateInitialRotPlot;

%% final periodic orbit

finalState = State{2}.finalConstraint(1:6);
finalPeriod = State{2}.period;
[stateFinalND, timeFinal] = periodicOrbit2Ephemeris(System{1}, ...
	finalState, finalPeriod, Option, EphemOption{2});

nSegment = length(timeFinal)-1;
stateFinalRotPlot = [];
JDFix.final = EphemOption{2}.JD;
finalnRev = EphemOption{2}.nRev;
finals = EphemOption{2}.s;
finaliFix = floor(finalnRev/2)*finals + 1; % index of the state to be fixed along the stack

for i = 1:nSegment
	ts = [timeFinal(i), timeFinal(i+1)];
	Y0 = stateFinalND(i, :);
	
	[T, Y] = ode113(@(t,y) ephemeris(t,y,JDFix.final, FrameSystem, Body), ts, Y0, opts);
	YRot = inert2rot([Y(:,1:3)*lstar, Y(:,4:6)*lstar/tstar], T*tstar, JDFix.final, FrameSystem);
	stateFinalRotPlot = [stateFinalRotPlot; YRot]; %#ok<AGROW>
	if i == finaliFix
		stateFix.final = Y0';
		figure(Option.plotNo(2))
		hold on
		axis equal
		grid on
		plot3(YRot(1,1), YRot(1,2), YRot(1,3), 'vr', 'linewidth', 2.0)
	end
end

figure(Option.plotNo(2))
plot3(stateFinalRotPlot(:,1), stateFinalRotPlot(:,2), stateFinalRotPlot(:,3), 'k')
stateRotPlot.final = stateFinalRotPlot;

end
