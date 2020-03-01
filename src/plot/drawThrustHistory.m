function drawThrustHistory(figureNo, System, State, Spacecraft)

nPhase = length(State);
figure(figureNo);
hold on

% isSame = isequaln(System{1}.dynamics, System{2}.dynamics);
isSame = 1;
if isSame
	for iPhase = 1:nPhase

		tstar = System{iPhase}.parameter.tstar;
		t = State{iPhase}.timeSegment;
		[thrustMaxND, ~, ~, thrustMaxD] = getSpacecraftInfo(Spacecraft{iPhase});
		controlMat = State{iPhase}.control;
		tday = t*tstar/60/60/24;
		Ttotal = controlMat(:,1)*thrustMaxD/thrustMaxND;
		Tx = Ttotal.*controlMat(:,2);
		Ty = Ttotal.*controlMat(:,3);
		Tz = Ttotal.*controlMat(:,4);

		h1 = stairs(tday, [Ttotal; Ttotal(end)], 'k', 'linewidth', 1.5);
		h2 = stairs(tday, [Tx; Tx(end)], 'r', 'linewidth', 1.5);
		h3 = stairs(tday, [Ty; Ty(end)], 'g', 'linewidth', 1.5);
		h4 = stairs(tday, [Tz; Tz(end)], 'b', 'linewidth', 1.5);
		
	end
else
	for iPhase = 1:nPhase

		tstar = System{iPhase}.parameter.tstar;
		t = State{iPhase}.timeSegment;
		[thrustMaxND, ~, ~, thrustMaxD] = getSpacecraftInfo(Spacecraft{iPhase});
		[~, ~, controlMat] = getStateControlMat(State{iPhase});
		if iPhase == 1
			toffset = 0
		else
			toffset = State{1}.timeSegment(end)*System{1}.parameter.tstar/60/60/24
		end
		tday = t*tstar/60/60/24 + toffset;
		Ttotal = controlMat(:,1)*thrustMaxD/thrustMaxND;
		Tx = Ttotal.*cos(controlMat(:,2)).*cos(controlMat(:,3));
		Ty = Ttotal.*sin(controlMat(:,2)).*cos(controlMat(:,3));
		Tz = Ttotal.*sin(controlMat(:,3));

		h1 = stairs(tday, [Ttotal; Ttotal(end)], 'k', 'linewidth', 1.5);
		h2 = stairs(tday, [Tx; Tx(end)], 'r', 'linewidth', 1.5);
		h3 = stairs(tday, [Ty; Ty(end)], 'g', 'linewidth', 1.5);
		h4 = stairs(tday, [Tz; Tz(end)], 'b', 'linewidth', 1.5);	
	
	end
end

legend([h1, h2, h3, h4], {'Tmag', 'Tx', 'Ty', 'Tz'});

end
