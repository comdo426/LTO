function InitialGuess = setInitialGuess(System, Setup)

% if Setup.transferType ... : If you are using manifolds. Room for improvements.

tau = [-1; -sqrt(495 +66*sqrt(15))/33; -sqrt(495 -66*sqrt(15))/33; 0; ...
	sqrt(495 -66*sqrt(15))/33; sqrt(495 +66*sqrt(15))/33; 1];
times = 0; % initialization for distinction from the function "times"
nodes = 0;
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-17);

nPhase = length(System);

InitialGuess = cell(nPhase,1);

for iPhase = 1:nPhase
	s = Setup.Phase{iPhase,1}.s;
	nRev = Setup.Phase{iPhase,1}.nRev;
	load(Setup.Phase{iPhase,1}.orbitName)
	isJC = strcmp(Setup.Phase{iPhase,1}.orbitSelect{1}, 'JC');
	nOrbit = length(times);
	mu = System{iPhase}.parameter.mu;
	if isJC
		targetJC = Setup.Phase{iPhase,1}.orbitSelect{2};
		targetJCArray = ones(1,nOrbit)*targetJC;
		JCArray = cell2mat(JC);
		diffJCArray = abs(targetJCArray - JCArray);
		[diffJC, orbitIndex] = min(diffJCArray);
		fprintf('Phase %d, Set JC: %d, Orbit JC: %d, difference: %d\n', ...
			iPhase, targetJC, JCArray(orbitIndex), diffJC);
		InitialGuess{iPhase, 1}.stateInitial = nodes{orbitIndex}(1,:); % initial condition
		period = times{orbitIndex}(end); % period
	else
		if iPhase == 1
			% 			InitialGuess{iPhase, 1}.stateInitial = [0.195328385711603; 0.0000000000000; 0; 0; 2.717614643431193 ; 0]'; % initial condition
			% 			period = 6.310163395436344; % period
			InitialGuess{iPhase, 1}.stateInitial = [1.115001963068551;                   0;   0.190879918445258;                   0;  -0.223386275592327;                   0]';
			period = 2*1.418376500453653;
			% 		fprintf('Wrong input for orbitSelect, under development');
		else
			% 			InitialGuess{iPhase, 1}.stateInitial = ...
			% 				[0.928542597235188; 0.866025403784438; 0; 0.110773126272870; -0.983802039797893; 0]'; % initial condition
			% 			period = 2*3.135521047051145; % period
			InitialGuess{iPhase, 1}.stateInitial = [0.840645749967222;                   0;   0.159210178951199;                   0;   0.261992933019641;                   0]';
			period = 2*1.349018321829528;
			
		end % Using JC if loop
	end
	
	
	
	if iPhase == 1
		tSegment = linspace(0, period, s+1)';
	else
		tSegment = timeFromPhase1 + linspace(0, period, s+1)';
	end
	dtSegment = tSegment(2) - tSegment(1);
	stateOneRev(1, :) = [InitialGuess{iPhase, 1}.stateInitial, 1];
	
	for iSegment = 1:s
		tNode = zeros(7,1);
		for itau = 1:7
			tNode(itau) = tSegment(iSegment) + dtSegment*(tau(itau)+1)/2;
		end % itau for loop
		[~, state] = ode113(@(t,y) CR3BP(t,y,mu,1,0), [tNode(1), tNode(3)], ...
			stateOneRev(3*(iSegment-1)+1, 1:6), opts);
		stateOneRev(3*(iSegment-1)+2, :) = [state(end, :), 1];
		[~, state] = ode113(@(t,y) CR3BP(t,y,mu,1,0), [tNode(3), tNode(5)], ...
			stateOneRev(3*(iSegment-1)+2, 1:6), opts);
		stateOneRev(3*(iSegment-1)+3, :) = [state(end, :), 1];
		if iSegment < s
			[~, state] = ode113(@(t,y) CR3BP(t,y,mu,1,0), [tNode(1), tNode(3)], ...
				stateOneRev(3*(iSegment-1)+3, 1:6), opts);
			stateOneRev(3*(iSegment-1)+4, :) = [state(end, :), 1];
		else
		end
	end % iSegment for loop
	
	stateMatrix = zeros(3*s*nRev + 1, 7);
	for iRev = 1:nRev
		if iRev < nRev
			stateMatrix(3*s*(iRev-1)+1:3*s*iRev, :) = stateOneRev;
		else
			stateMatrix(3*s*(iRev-1)+1:3*s*iRev+1, :) = ...
				[stateOneRev; stateOneRev(1, :)];
		end % if loop
	end % iRev for loop
	stateArray = reshape(stateMatrix', [7*(3*s*nRev + 1), 1]);
	
	time = zeros(s*nRev+1, 1);
	for iRev = 1:nRev
		if iRev < nRev
			time(s*(iRev-1)+1: s*iRev) = tSegment(1:end-1)+(iRev-1)*period;
		else
			time(s*(iRev-1)+1: s*iRev+1) = ...
				[tSegment(1:end-1)+(iRev-1)*period; tSegment(1) + nRev*period];
			timeFromPhase1 = time(end);
		end
	end
	
	controlMatrix = repmat([1e-10, 0, 0], [s*nRev 1]); % play with numbers if you want
	controlArray = reshape(controlMatrix', [3*s*nRev, 1]);
	
	% Converted to column vector, transpose needed to make the dimensions right
	InitialGuess{iPhase, 1}.state = stateArray;
	InitialGuess{iPhase, 1}.control = controlArray;
	InitialGuess{iPhase, 1}.timeSegment = time;
	InitialGuess{iPhase, 1}.nSegment = s*nRev; % Remember that segment requires two nodes
	InitialGuess{iPhase, 1}.t0 = 0; % Time offset of the rotating frame
	
	if iPhase == 1
		InitialGuess{iPhase, 1}.initialConstraint = ...
			[InitialGuess{iPhase, 1}.stateInitial, 1]';
	else
		if iPhase == nPhase
			InitialGuess{iPhase, 1}.finalConstraint = ...
				stateOneRev(1, 1:6)'; % open to change, only constrain 6 states
		end
	end
	
	InitialGuess{iPhase, 1}.slack = [];
	
	if Setup.plot{1}
		stateSegmentMatrix = zeros(s*nRev+1,7);
		for iSegment = 1:InitialGuess{iPhase, 1}.nSegment+1
			stateSegmentMatrix(iSegment, :) = stateMatrix(3*(iSegment-1)+1, :);
		end
		
		figure(Setup.plot{2})
		commonAxisSetting;
		CR3BPAxisSetting;
		plot3(stateSegmentMatrix(:,1), stateSegmentMatrix(:,2), stateSegmentMatrix(:,3), 'k.', 'linewidth', 2.0)
		plot3(stateMatrix(:,1), stateMatrix(:,2), stateMatrix(:,3), 'b-', 'linewidth', 1.0);
	end
	
	clear stateOneRev;
end % iPhase for loop

if Setup.plot{1}
	figure(Setup.plot{2})
	earthPlot = drawEarth(mu);
	moonPlot = drawMoon(mu);
	lpPlot = drawLagrangianPoints(mu, Setup.plot{3});
end
[initialPlot, ~, finalPlot] = drawEndPoints(Setup.plot{2}, InitialGuess);

figure(Setup.plot{2})
legend([initialPlot, finalPlot, ...
	moonPlot, lpPlot], {'Initial', 'Final', ...
	 'Moon', 'L_i'}, 'fontsize', 12);




end