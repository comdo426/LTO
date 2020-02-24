function StateEphemeris = ...
	getEphemerisTransfer(System, State, Spacecraft, ...
	Option, EphemOption, JDFix, stateFix)
%GETEPHEMERISTRANSFER - get transfer trajectory in ephemeris
%
%  Syntax:
%     StateEphemeris = ...
%		GETEPHEMERISTRANSFER(System, State, Spacecraft, ...
%		Option, EphemOption, JDFix, stateFix)
%
%  Description:
%     Stack periodic orbits by nRev / s segments and converge them in ephemeris
%     to get the boundary constraints in ephemeris
%
%  Outputs:
%		StateEphemeris - structure that contains trajectory converged in
%		ephemeris.
%			stateInertMat - state in ephemeris, inertial, non-dimensional unit
%			stateConvergedRotPlot - just for plot, rotational state in
%			non-dimensional unit
%			timeConverged - time converged in non-dimensional units
%			controlMat - control information			
%
%  See also: TRANSFER2EPHEMERISMAIN, GETEPHEMERISBOUNDARYCONST
%
%   Author: Beom Park
%   Date: 23-Feb-2020; Last revision: 23-Feb-2020

% merge the states 
[~, stateSegment1Mat, control1Mat] = getStateControlMat(State{1});
[~, stateSegment2Mat, control2Mat] = getStateControlMat(State{2});
stateSegmentMerged = [stateSegment1Mat; stateSegment2Mat(2:end, :)];
nSegment = State{1}.nSegment + State{2}.nSegment;
tJump = State{1}.timeSegment(end);
timeMerged = [State{1}.timeSegment; tJump + State{2}.timeSegment(2:end)];
control2Mat(:,2) = control2Mat(:,2) + tJump;
controlMerged = [control1Mat; control2Mat];

% get the initial phasing angle 
secPastJ2000 = (JDFix.initial-cspice_j2000)*60*60*24;
moonState = cspice_spkezr('MOON', secPastJ2000, Option.FrameSystem.frame, ...
	'NONE', 'EARTH');
moonPos = moonState(1:3);
xUnit = [1, 0, 0];
theta0 = acos(dot(moonPos, xUnit)/norm(moonPos)/norm(xUnit));
if moonPos(2) < 0
	theta0 = -theta0;
end
% Add the initial phasing angle
controlMerged(:,2) = controlMerged(:,2) + theta0;

% FrameSystem, Body, opts, getSpacecraftInfo
FrameSystem = Option.FrameSystem;
tstar = FrameSystem.tstar;
lstar = FrameSystem.lstar;
Body = Option.Body;
opts = Option.integrate;
[TmaxND, IspND, g0ND, ~, ~] = getSpacecraftInfo(Spacecraft{1});
mass = stateSegmentMerged(:,7);
stateRotPosVel = stateSegmentMerged(:, 1:6);
stateInertPosVel = rot2inert(stateRotPosVel, timeMerged*tstar, JDFix.initial, ...
	FrameSystem);
stateInertPosVelND = ...
	[stateInertPosVel(:,1:3)/lstar, stateInertPosVel(:,4:6)/lstar*tstar];
stateInert = [stateInertPosVelND, mass];

% Propagate the initial trajectory in ephemeris to see the initial deviation
stateTransferRotPlot = [];
for i = 1:nSegment
	Y0 = stateInert(i, :);
	u = controlMerged(i,:);
	ts = [timeMerged(i), timeMerged(i+1)];
	[T, Y] = ode113(@(t,y)...
		ephemerisLT(t, y, u, JDFix.initial, IspND, g0ND, FrameSystem, Body), ...
		ts, Y0, opts);	
	YRot = inert2rot([Y(:,1:3)*lstar, Y(:,4:6)*lstar/tstar], T*tstar, ...
		JDFix.initial, FrameSystem);
	stateTransferRotPlot = [stateTransferRotPlot; [YRot, Y(:,7)]]; %#ok<AGROW>
end
figure(Option.plotNo(3))
axis equal
grid on
hold on
plot3(stateTransferRotPlot(:,1), stateTransferRotPlot(:,2), ... 
	stateTransferRotPlot(:,3), 'c', 'linewidth', 1.0)

% Call the solver
Problem = setProblemEphemerisLT(stateInert, timeMerged, controlMerged, ...
	JDFix, stateFix, IspND, g0ND, TmaxND, FrameSystem, Body, Option);
converged = newtonRaphson(Problem);

% Parse the converged solution to draw plots
nState = 7*nSegment+7;
nTime = nSegment+1;
stateConvergedVec = converged(1:nState);
stateInertConvergedND = reshape(stateConvergedVec, [7, nTime])'; ...
	%Transpose to match the dimensions
timeConverged = converged(nState+nSegment+1:nState+nSegment+nTime);

controlConvergedVec = ...
	converged(nState+nSegment+nTime+1:nState+nSegment+nTime+3*nSegment);
controlMat = reshape(controlConvergedVec, [3, nSegment])';

stateConvergedRotPlot = [];

for i = 1:nSegment
	Y0 = stateInertConvergedND(i, :);
	u = controlMat(i,:);
	ts = [timeConverged(i), timeConverged(i+1)];
	[T, Y] = ode113(@(t,y) ...
		ephemerisLT(t, y, u, JDFix.initial, IspND, g0ND, FrameSystem, Body), ...
		ts, Y0, opts);
	
	YRot = inert2rot([Y(:,1:3)*lstar, Y(:,4:6)*lstar/tstar], ...
		T*tstar, JDFix.initial, FrameSystem);
	stateConvergedRotPlot = [stateConvergedRotPlot; [YRot, Y(:,7)]]; %#ok<AGROW>
end

% Define outputs
StateEphemeris.stateInertMat = stateInertConvergedND;
StateEphemeris.stateConvergedRotPlot = stateConvergedRotPlot;
StateEphemeris.timeConverged = timeConverged;
StateEphemeris.controlMat = controlMat;

% Draw plots

figure(Option.plotNo(4))
hold on
grid on
axis equal
plot3(stateConvergedRotPlot(:,1), stateConvergedRotPlot(:,2), ...
	stateConvergedRotPlot(:,3), 'c', 'linewidth', 2.0)

tday = timeConverged*tstar/60/60/24;
Ttotal = controlMat(:,1)/TmaxND*Spacecraft{1}.thrustMaxD;
Tx = Ttotal.*cos(controlMat(:,2)).*cos(controlMat(:,3));
Ty = Ttotal.*sin(controlMat(:,2)).*cos(controlMat(:,3));
Tz = Ttotal.*sin(controlMat(:,3));

figure(Option.plotNo(5))
hold on
grid on
h1 = stairs(tday, [Ttotal; Ttotal(end)], 'k', 'linewidth', 1.5);
h2 = stairs(tday, [Tx; Tx(end)], 'r', 'linewidth', 1.5);
h3 = stairs(tday, [Ty; Ty(end)], 'g', 'linewidth', 1.5);
h4 = stairs(tday, [Tz; Tz(end)], 'b', 'linewidth', 1.5);
legend([h1, h2, h3, h4], {'T', 'T_x', 'T_y', 'T_z'})
xlabel('time(days)')
ylabel('thrust(N)')

end