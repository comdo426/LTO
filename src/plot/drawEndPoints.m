function [initial, interPhase, final] = drawEndPoints(figureNo, System, State)



figure(figureNo)
commonAxisSetting


nPhase = length(State);


% isSame = isequaln(System{1}.dynamics, System{2}.dynamics);
isSame = 1;

if isSame
	
	%% initialpoint
	
	initialPoint = State{1}.initialConstraint;
	initial = plot3(initialPoint(1), initialPoint(2), initialPoint(3), '^k', 'linewidth', 2.0);
	
	%% finalpoint
	
	finalPoint = State{nPhase}.finalConstraint;
	final = plot3(finalPoint(1), finalPoint(2), finalPoint(3), 'vk', 'linewidth', 2.0);
	
	%% intermediatePoint
	
	if nPhase > 1
		interPoint = nan(nPhase-1, 3);
		
		
		for iPhase = 1:nPhase-1
			interPoint(iPhase, :) = State{iPhase+1}.state(1, 1:3);
		end
		
		
		interPhase = plot3(interPoint(:, 1), interPoint(:, 2), interPoint(:, 3), 'dk', 'linewidth', 2.0);
	else
		interPhase = [];
	end
	
else
	
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
	JD0 = State{1}.JD0;
	figure(figureNo+100)
	commonAxisSetting;
	
	%% initialpoint
	
	figure(figureNo)
	initialPoint = State{1}.state(1:3);
	plot3(initialPoint(1), initialPoint(2), initialPoint(3), '^k', 'linewidth', 2.0);
	
	secPastJ2000 = (JD0-cspice_j2000)*60*60*24;
	initialInertEarth = rot2inert(State{1}.state(1:6)', 0, JD0, EM);
	EarthInitial = cspice_spkezr('EARTH', secPastJ2000, 'J2000', 'NONE', 'SUN');
	
	initialInert = initialInertEarth + EarthInitial';
	initialInertND = [initialInert(1:3)/lstarS, initialInert(4:6)/lstarS*tstarS];
	
	figure(figureNo+100)
	initial = plot3(initialInertND(1), initialInertND(2), initialInertND(3), '^k', 'linewidth', 2.0);
	
	
	%% finalpoint
	
	finalPoint = State{2}.state(end-6:end-4);
	figure(figureNo+100)
	final = plot3(finalPoint(1), finalPoint(2), finalPoint(3), 'vk', 'linewidth', 2.0);
	
	%% intermediatepoint
	
	interPoint = State{2}.state(1:6);
	figure(figureNo+100)
	interPhase = plot3(interPoint(1), interPoint(2), interPoint(3), 'dk', 'linewidth', 2.0);
	
	secPastJ2000 = State{1}.timeSegment(end)*tstar + (JD0-cspice_j2000)*60*60*24;
	EarthInter = cspice_spkezr('EARTH', secPastJ2000, 'J2000', 'NONE', 'SUN');
	interPointD = [interPoint(1:3)'*lstarS, interPoint(4:6)'*lstarS/tstarS];
	interInertEarth = interPointD - EarthInter';
	interRot = inert2rot(interInertEarth, State{1}.timeSegment(end)*tstar, JD0, ...
		EM);
	
	figure(figureNo)
	plot3(interRot(1), interRot(2), interRot(3), 'dk', 'linewidth', 2.0);
	
end
end