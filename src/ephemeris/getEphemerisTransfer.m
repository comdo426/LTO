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
%   Date: 03-Mar-2020; Last revision: 03-Mar-2020

% merge the states
state1Mat = State{1}.state;
[~, idx1] = getIndexSegment(State{1});
stateSegment1Mat = state1Mat(idx1, :);
control1Mat = State{1}.control;

state2Mat = State{2}.state;
[~, idx2] = getIndexSegment(State{2});
stateSegment2Mat = state2Mat(idx2, :);
control2Mat = State{2}.control;

stateSegmentMerged = [stateSegment1Mat; stateSegment2Mat; state2Mat(end, :)];
nSegment = State{1}.nSegment + State{2}.nSegment;
% From the latest version, tJump is internally taken care of
tJump = 0;
timeMerged = [State{1}.timeSegment; tJump + State{2}.timeSegment(2:end)];
% control2Mat(:,2) = control2Mat(:,2) + tJump;
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
% Incorporate the initial phasing angle to the hey
C = [cos(theta0), -sin(theta0);
	sin(theta0), cos(theta0)];
for i = 1:nSegment
	controlMerged(i,2:3) = (C*controlMerged(i,2:3)')';
end

% FrameSystem, Body, opts, getSpacecraftInfo
FrameSystem = Option.FrameSystem;
mu = FrameSystem.mu;
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
stateTransferCR3BPPlot = [];
for i = 1:nSegment
	Y0 = stateInert(i, :);
	u = controlMerged(i,:);
	ts = [timeMerged(i), timeMerged(i+1)];
	[T, Y] = ode113(@(t,y)...
		ephemerisLT(t, y, u, JDFix.initial, IspND, g0ND, FrameSystem, Body), ...
		ts, Y0, opts);
	YRot = inert2rot([Y(:,1:3)*lstar, Y(:,4:6)*lstar/tstar], T*tstar, ...
		JDFix.initial, FrameSystem);
	[T2, Y2] = ode113(@(t,y) CR3BPLT(t,y,u,mu, IspND, g0ND, 0), ts, stateSegmentMerged(i,:), opts);
	stateTransferRotPlot = [stateTransferRotPlot; [YRot, Y(:,7)]]; %#ok<AGROW>
	stateTransferCR3BPPlot = [stateTransferCR3BPPlot; Y2];  %#ok<AGROW>
end
figure(Option.plotNo(3))
axis equal
grid on
hold on
plot3(stateTransferRotPlot(:,1), stateTransferRotPlot(:,2), ...
	stateTransferRotPlot(:,3), 'r', 'linewidth', 1.0)
plot3(stateTransferCR3BPPlot(:,1), stateTransferCR3BPPlot(:,2), ...
	stateTransferCR3BPPlot(:,3), 'k', 'linewidth', 1.0)
% 
% save('TESTInsideGetEphem2.mat')
% error(' ');

%% Call the solver(Multiple Shooter version)

isMultipleShooting = Option.isMultipleShooting;
if isMultipleShooting
	
	Option.newton.maxIteration = 30;
	
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
		converged(nState+nSegment+nTime+1:nState+nSegment+nTime+4*nSegment);
	controlMat = reshape(controlConvergedVec, [4, nSegment])';
	
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
		stateConvergedRotPlot(:,3), 'r', 'linewidth', 0.5)
	
	tday = timeConverged*tstar/60/60/24;
	Ttotal = controlMat(:,1)/TmaxND*Spacecraft{1}.thrustMaxD;
	Tx = Ttotal.*controlMat(:,2);
	Ty = Ttotal.*controlMat(:,3);
	Tz = Ttotal.*controlMat(:,4);
	
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
	
else
	
	%% Call the solver(Collocation version)
	
	Option.JDFix = JDFix;
	Option.stateFix = stateFix;
	Option.newton.maxIteration = 50;
	Option.C = C;
	
	Collocation = setCollocation;
	
	[Problem, State] = setProblemEphemerisLTCollocation(System, State, Spacecraft, Option, ...
		Collocation);
	
	feasibleWithSlack = newtonRaphson(Problem);
	feasibleVec = deleteSlackVariable(feasibleWithSlack, State, Option);
	State = updateState(State, feasibleVec);
	
	% Phase 1
	
	stateTransferRotPlotCollocation = [];
	
	nPhase = length(System);
	for iPhase = 1:nPhase
		stateMat = State{iPhase}.state;
		[~, indexSegment] = getIndexSegment(State{iPhase});
		stateSegmentMat = stateMat(indexSegment, :);
		controlMat = State{iPhase}.control;
		nSegment = State{iPhase}.nSegment;
		t = State{iPhase}.timeSegment;
		tVar = State{iPhase}.timeVariable;
		
		%
		for i = 1:nSegment
			Y0 = stateSegmentMat(i, :);
			u = controlMat(i,:);
			ts = [t(i), t(i+1)];
			[T, Y] = ode113(@(t,y)...
				ephemerisLT(t, y, u, JDFix.initial, IspND, g0ND, FrameSystem, Body), ...
				ts, Y0, opts);
			YRot = inert2rot([Y(:,1:3)*lstar, Y(:,4:6)*lstar/tstar], T*tstar, ...
				JDFix.initial, FrameSystem);
			stateTransferRotPlotCollocation = [stateTransferRotPlotCollocation; ...
				[YRot, Y(:,7)]]; %#ok<AGROW>
			if i < nSegment
				error = Y(end, :) - stateSegmentMat(i+1,:);
				errorMat(iPhase, i, :) = error;
				errorNorm(iPhase,i) = norm(error);
			end
		end
		
		stateRot = inert2rot([stateMat(:,1:3)*lstar, stateMat(:,4:6)*lstar/tstar], ...
			tVar*tstar, JDFix.initial, FrameSystem);
		
% 		figure(Option.plotNo(4))
% 		axis equal
% 		hold on
% 		plot3(stateRot(:,1), stateRot(:,2), stateRot(:,3), 'r', 'linewidth', 1.5)
% 		StateEphemeris.State = State;
% 		StateEphemeris.stateConvergedRotPlot = stateTransferRotPlotCollocation;
% 		
		
	end
	%
	% figure(Option.plotNo(4))
	% axis equalload
	% grid on
	% hold on
	% plot3(stateTransferRotPlotCollocation(:,1), stateTransferRotPlotCollocation(:,2), ...
	% 	stateTransferRotPlotCollocation(:,3), 'r', 'linewidth', 1.0)
	isOptimize = LToptget(Option.LTO, 'Optimizer');
	
	if isOptimize
		
		cprintf(-[1, 0, 0], 'Optimize: 1\n');
		[x0, funcs, options, State] = setProblemipoptEphemeris(System, State, Spacecraft, Option, Collocation);
		[optimizedWithSlack, ~] = ipopt_auxdata(x0, funcs, options);
		optimizedVec = deleteSlackVariable(optimizedWithSlack, State, Option);
		State = updateState(State, optimizedVec);
		stateTransferRotPlotCollocationOptimized = [];
		for iPhase = 1:nPhase
			stateMat = State{iPhase}.state;
			[~, indexSegment] = getIndexSegment(State{iPhase});
			stateSegmentMat = stateMat(indexSegment, :);
			controlMat = State{iPhase}.control;
			nSegment = State{iPhase}.nSegment;
			t = State{iPhase}.timeSegment;
			tVar = State{iPhase}.timeVariable;
			
			%
			for i = 1:nSegment
				Y0 = stateSegmentMat(i, :);
				u = controlMat(i,:);
				ts = [t(i), t(i+1)];
				[T, Y] = ode113(@(t,y)...
					ephemerisLT(t, y, u, JDFix.initial, IspND, g0ND, FrameSystem, Body), ...
					ts, Y0, opts);
				YRot = inert2rot([Y(:,1:3)*lstar, Y(:,4:6)*lstar/tstar], T*tstar, ...
					JDFix.initial, FrameSystem);
				stateTransferRotPlotCollocationOptimized = [stateTransferRotPlotCollocationOptimized; ...
					[YRot, Y(:,7)]]; %#ok<AGROW>
				if i < nSegment
					error = Y(end, :) - stateSegmentMat(i+1,:);
					errorOptMat(iPhase, i, :) = error;
					errorOptNorm(iPhase,i) = norm(error);
				end
					
			end
			
			stateRot = inert2rot([stateMat(:,1:3)*lstar, stateMat(:,4:6)*lstar/tstar], ...
				tVar*tstar, JDFix.initial, FrameSystem);
			
			figure(Option.plotNo(4))
			axis equal
			hold on
			plot3(stateRot(:,1), stateRot(:,2), stateRot(:,3), 'r', 'linewidth', 1.5)
			StateEphemeris.State = State;
			StateEphemeris.stateConvergedRotPlot = stateTransferRotPlotCollocation;
			StateEphemeris.stateConvergedRotPlotOptimized = stateTransferRotPlotCollocationOptimized;
			
			
		end
	end
	
	save('TEST')
	
end