nState = 7*nSegment+7;
nTime = nSegment+1;
stateConvergedVec = converged(1:nState);
stateInertConvergedND = reshape(stateConvergedVec, [7, nTime])'; ...
	%Transpose to match the dimensions
timeConverged = converged(nState+nSegment+1:nState+nSegment+nTime);
timeConvergedSecond = timeConverged*tstar;

controlConvergedVec = converged(nState+nSegment+nTime+1:nState+nSegment+nTime+3*nSegment);
controlMat = reshape(controlConvergedVec, [3, nSegment])';
% 
% stateInertConverged = [stateInertConvergedND(:,1:3)*lstar, ...
% 	stateInertConvergedND(:,4:6)*lstar/tstar];

stateInitialRotPlot = [];

for i = 1:nSegment
	Y0 = stateInertConvergedND(i, :);
	u = controlMat(i,:);
	ts = [timeConverged(i), timeConverged(i+1)];
	[T, Y] = ode113(@(t,y) ephemerisLT(t, y, u, initialJDFix, IspND, g0ND, EM, Body), ...
		ts, Y0, opts);
	
	YRot = inert2rot([Y(:,1:3)*lstar, Y(:,4:6)*lstar/tstar], T*tstar, initialJDFix, EM);
	stateInitialRotPlot = [stateInitialRotPlot; [YRot, Y(:,7)]];
end

figure(5)
hold on
grid on
axis equal
plot3(stateInitialRotPlot(:,1), stateInitialRotPlot(:,2), stateInitialRotPlot(:,3), 'c', 'linewidth', 2.0)

tday = timeConverged*tstar/60/60/24;
Ttotal = controlMat(:,1)/TmaxND*Spacecraft{1}.thrustMaxD;
Tx = Ttotal.*cos(controlMat(:,2)).*cos(controlMat(:,3));
Ty = Ttotal.*sin(controlMat(:,2)).*cos(controlMat(:,3));
Tz = Ttotal.*sin(controlMat(:,3));

figure(4)
hold on
h1 = stairs(tday, [Ttotal; Ttotal(end)], 'k', 'linewidth', 1.5);
h2 = stairs(tday, [Tx; Tx(end)], 'r', 'linewidth', 1.5);
h3 = stairs(tday, [Ty; Ty(end)], 'g', 'linewidth', 1.5);
h4 = stairs(tday, [Tz; Tz(end)], 'b', 'linewidth', 1.5);