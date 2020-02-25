function [earthPlot, moonPlot, lpPlot] = ...
	drawThrustArc(figureNo, LPointVec, System, State, Spacecraft)
%DRAWTHRUSTARC - draw the trajectory arc after convergence
%
%  Syntax:
%     [earthPlot, moonPlot, lpPlot] = ...
%		DRAWTHRUSTARC(figureNo, LPointVec, System, State, Spacecraft)
%
%  Description:
%     draws the trajectory arc after convergence. It currently supports EM CR3BP
%     or EM-CR3BP + Sun centered conic. It needs to be generalized for other
%     dynamics model. 
%
%  Inputs:
%		figureNo - number of the figure to draw on
%		LPointVec - vector that contains the number of Lagrangain points
%
%  Outputs:
%     earthPlot - plot handle for Earth
%		moonPlot - plot handle for moon
%		lpPlot - ploot handle for Lagrangian points
%
%	See also: DRAWENDPOINTS, LTOMAIN
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 17-Feb-2020

% mu = System{1}.parameter.mu;
nPhase = length(State);

% First figure for the rotating frame
figure(figureNo)
commonAxisSetting;
CR3BPAxisSetting;

% Second figure for the inertial frame
figure(figureNo+100)
commonAxisSetting;

earthPlot = [];
moonPlot = [];
lpPlot = [];
% earthPlot = drawEarth(mu);
% moonPlot = drawMoon(mu);
% lpPlot = drawLagrangianPoints(mu, LPointVec);

% NOTE: temporary change to make it work for CR3BP - 2BP
% isSame = isequaln(System{1}.dynamics, System{2}.dynamics);
isSame = 1;

if isSame
	
	for iPhase = 1:nPhase
		[~, nSegment, ~, ~, ~, ~, ~] = getPhaseStateInfo(State{iPhase});
		stateMat = State{iPhase}.state;
		controlMat = State{iPhase}.control;
		indexSegment = getIndexSegment(State{iPhase});
		stateSegmentMat = stateMat(indexSegment, :);
		for i = 1:nSegment
			if controlMat(i, 1) > 0.1*Spacecraft{iPhase}.thrustMaxND
				color = [1 0 0];
			else
				color = [0 0 1];
			end
			figure(figureNo)
			plot3(stateMat(4*(i-1)+1:4*i,1), stateMat(4*(i-1)+1:4*i,2), ...
				stateMat(4*(i-1)+1:4*i,3), 'Color', color, ...
				'linewidth', 1.5);
		end
		figure(figureNo)
		plot3(stateSegmentMat(:,1), stateSegmentMat(:,2), stateSegmentMat(:,3), ...
			'k.', 'linewidth', 1.0);
	end
	
else % only works for two phases
	
	%% first phase: rotational frame
	
	EM.P1 = 'EARTH';
	EM.P2 = 'MOON';
	EM.mu1 = System{1}.parameter.mu1;
	EM.mu2 = System{1}.parameter.mu2;
	EM.centralBody = 'EARTH';
	EM.frame = 'J2000';
	mu = System{1}.parameter.mu;
	lstar = System{1}.parameter.lstar;
	tstar = System{1}.parameter.tstar;
	lstarS = System{2}.parameter.lstar;
	tstarS = System{2}.parameter.tstar;
	[t, nSegment, ~, ~, ~, ~, ~] = getPhaseStateInfo(State{1});
	[stateMat, stateSegmentMat, controlMat] = ...
		getStateControlMat(State{1});
	
	JD0 = State{1}.JD0;
	
	% rotation
	
	t = t-t(1);
	JDf = JD0 + t(end)*tstar/60/60/24;
	
	tau = [-1, -0.830223896278567, -0.468848793470714, 0, ...
		0.468848793470714, 0.830223896278567, 1];
	
	tAug = nan(3*nSegment+1,1);
	for i=1:nSegment
		dt = t(i+1)-t(i);
		tAug(3*(i-1)+1) = t(i);
		tAug(3*(i-1)+2) = t(i) + dt*(tau(3)+1)/2;
		tAug(3*(i-1)+3) = t(i) + dt*(tau(5)+1)/2;
	end
	tAug(3*nSegment+1) = t(nSegment+1);
	etVec = (tAug*tstar + (JD0-cspice_j2000)*60*60*24)';
	flybyEarthSun = cspice_spkezr('EARTH', etVec, 'J2000', 'NONE', 'SUN');
	
	flybyInertEarth = rot2inert(stateMat, tAug*tstar, JD0, EM);
	flybyInertSun = flybyInertEarth + flybyEarthSun';
	flybyInertSunND = [flybyInertSun(:,1:3)/lstarS, ...
		flybyInertSun(:,4:6)/lstarS*tstarS];
	
	etVecSeg = (t*tstar + (JD0-cspice_j2000)*60*60*24)';
	flybyEarthSunSeg = cspice_spkezr('EARTH', etVecSeg, 'J2000', 'NONE', 'SUN');
	
	flybyInertEarthSeg = rot2inert(stateSegmentMat, t*tstar, JD0, EM);
	flybyInertSunSeg = flybyInertEarthSeg + flybyEarthSunSeg';
	flybyInertSunSegND = [flybyInertSunSeg(:, 1:3)/lstarS, ...
		flybyInertSunSeg(:,4:6)/lstarS*tstarS];
	
	
	
	for i = 1:nSegment
		if controlMat(i, 1) > 0.1*Spacecraft{1}.thrustMaxND
			color = [1 0 0];
		else
			color = [0 0 1];
		end
		figure(figureNo)
		plot3(stateMat(3*(i-1)+1:3*i+1,1), stateMat(3*(i-1)+1:3*i+1,2), ...
			stateMat(3*(i-1)+1:3*i+1,3), 'Color', color, ...
			'linewidth', 1.5);
		
		figure(figureNo+1)
		hold on
		axis equal
		plot3(stateMat(3*(i-1)+1:3*i+1,4), stateMat(3*(i-1)+1:3*i+1,5), ...
			stateMat(3*(i-1)+1:3*i+1,6), 'k', ...
			'linewidth', 1.5);
		figure(figureNo+100)
		% do the transformation here and move it to inertial frame
		plot3(flybyInertSunND(3*(i-1)+1:3*i+1,1), flybyInertSunND(3*(i-1)+1:3*i+1,2), ...
			flybyInertSunND(3*(i-1)+1:3*i+1,3), 'Color', color, ...
			'linewidth', 1.5);
	end
	
	% 	figure(figureNo)
	% 	plot3(stateSegmentMat(:,1), stateSegmentMat(:,2), stateSegmentMat(:,3), ...
	% 		'k.', 'linewidth', 1.0);
	% 	figure(figureNo+100)
	% 	plot3(flybyInertSunSegND(:,1), flybyInertSunSegND(:,2), flybyInertSunSegND(:,3), ...
	% 		'k.', 'linewidth', 1.0);
	%
	
	
	%% second phase: inertial frame
	
	[t2, nSegment, ~, ~, ~, ~, ~] = getPhaseStateInfo(State{2});
	[stateMat, stateSegmentMat, controlMat] = ...
		getStateControlMat(State{2});
	
	tAug2 = nan(3*nSegment+1,1);
	for i=1:nSegment
		dt = t2(i+1)-t2(i);
		tAug2(3*(i-1)+1) = t2(i);
		tAug2(3*(i-1)+2) = t2(i) + dt*(tau(3)+1)/2;
		tAug2(3*(i-1)+3) = t2(i) + dt*(tau(5)+1)/2;
	end
	tAug2(end) = t2(end);
	
	etVec2 = (tAug2*tstarS + (JDf-cspice_j2000)*60*60*24)';
	etVecSeg2 = (t2*tstarS + (JDf - cspice_j2000)*60*60*24);
	
	legEarthSun = cspice_spkezr('EARTH', etVec2, 'J2000', 'NONE', 'SUN');
	stateMatD = [stateMat(:,1:3)*lstarS, stateMat(:,4:6)*lstarS/tstarS];
	legInertEarth = stateMatD - legEarthSun';
	legRotEarth = inert2rot(legInertEarth, tAug2*tstarS, JDf, EM);
	
	legEarthSunSeg = cspice_spkezr('EARTH', etVecSeg2, 'J2000', 'NONE', 'SUN');
	stateSegmentMatD = [stateSegmentMat(:,1:3)*lstarS, stateSegmentMat(:,4:6)*lstarS/tstarS];
	legInertEarthSeg = stateSegmentMatD - legEarthSunSeg';
	legRotEarthSeg = inert2rot(legInertEarthSeg, t2*tstarS, JDf, EM);
	
	for i = 1:nSegment
		if controlMat(i, 1) > 0.1*Spacecraft{1}.thrustMaxND
			color = [1 0 0];
		else
			color = [0 0 1];
		end
		figure(figureNo+100)
		plot3(stateMat(3*(i-1)+1:3*i+1,1), stateMat(3*(i-1)+1:3*i+1,2), ...
			stateMat(3*(i-1)+1:3*i+1,3), '--', 'Color', color, ...
			'linewidth', 1.5);
		
		% do the transformation here and move it to rotating frame
		figure(figureNo)
		plot3(legRotEarth(3*(i-1)+1:3*i+1,1), legRotEarth(3*(i-1)+1:3*i+1,2), ...
			legRotEarth(3*(i-1)+1:3*i+1,3), '--', 'Color', color, ...
			'linewidth', 1.5);
		figure(figureNo+1)
		plot3(legRotEarth(3*(i-1)+1:3*i+1,4), legRotEarth(3*(i-1)+1:3*i+1,5), ...
			legRotEarth(3*(i-1)+1:3*i+1,6), '--k', ...
			'linewidth', 1.5);
		
		
	end
	% 	figure(figureNo+100)
	% 	plot3(stateSegmentMat(:,1), stateSegmentMat(:,2), stateSegmentMat(:,3), ...
	% 		'k.', 'linewidth', 1.0);
	%
	% do the transformation here and move it to rotating frame
	% 	figure(figureNo)
	% 	plot3(legRotEarthSeg(:,1), legRotEarthSeg(:,2), legRotEarthSeg(:,3), ...
	% 		'k.', 'linewidth', 1.0);
	%
	%
	
	save('plotTest')
	
end
end