function InitialGuessPhase = setInitialGuessPhase(s, nRev, stateInitial, period, ...
	SystemPhase)
%propagateStateInitial - sets InitialGuess cell for each phase by propagation.

tOneRev = linspace(0, period, s+1)';
intRev = floor(nRev); % integer rev number of the given nRev
isInt = (intRev == nRev); % is Integer boolean

if nRev < 1
	error('At least one revolution needed for each phase')
end

stateOneRev = propagateStateInitial(s, tOneRev, stateInitial, false, SystemPhase);

if isInt % for integer revolutions
	nSegment = s*nRev;
else % for non-integer revolutions
	extraRev = nRev - intRev; % What is the 0.xxxx revolution remaining?
	sExtra = ceil(s*extraRev); % additional number of the segments for the extra rev
	tExtra = linspace(0, period*extraRev, sExtra+1)';
	stateExtra = propagateStateInitial(sExtra, tExtra, stateInitial, true, SystemPhase);
	nSegment = s*intRev + sExtra;
end

% state formulation
stateMatrix = nan(4*nSegment, 7);
for iRev = 1:intRev
	if iRev < intRev
		stateMatrix(4*s*(iRev-1)+1:4*s*iRev, :) = stateOneRev;
	elseif isInt
		stateMatrix(4*s*(iRev-1)+1:4*s*iRev, :) = stateOneRev;
	else
		stateMatrix(4*s*(iRev-1)+1:4*s*iRev, :) = stateOneRev;
		stateMatrix(4*s*iRev+1:4*nSegment, :) = stateExtra;
	end % if loop
end % iRev for loop
stateArray = reshape(stateMatrix', [7*4*nSegment, 1]);

% time formulation
time = nan(nSegment+1, 1);
for iRev = 1:intRev
	if iRev < intRev
		time(s*(iRev-1)+1: s*iRev) = tOneRev(1:end-1)+(iRev-1)*period;
	elseif isInt
		time(s*(iRev-1)+1: s*iRev+1) = ...
			[tOneRev(1:end-1)+(iRev-1)*period; intRev*period];
	else
		time(s*(iRev-1)+1: s*iRev) = tOneRev(1:end-1)+(iRev-1)*period;
		time(s*iRev+1: nSegment+1) = tExtra + intRev*period;		
	end
end

% time augmentation
timeVariable = nan(4*nSegment, 1);
timeDefect = nan(3*nSegment, 1);
Collocation = setCollocation;
for i = 1:nSegment
	dt = time(i+1) - time(i);
	timeAug = time(i) + dt*Collocation.tauRatio;
	timeVariable(4*(i-1)+1:4*i) = timeAug(1:2:7);
	timeDefect(3*(i-1)+1:3*i) = timeAug(2:2:6);
end


% control formulation

controlMatrix = nan(nSegment, 4); % thrust magnitude, thrust vector components
for i = 1:nSegment
	vel = stateMatrix(i, 4:6);
	controlMatrix(i, 1) = 1e-10; 
	controlMatrix(i, 2:4) = vel/norm(vel);
	C = [cos(time(i)), sin(time(i));
		-sin(time(i)), cos(time(i))];
	controlMaxtrix(i, 2:3) = (C*controlMatrix(i, 2:3)')';
end
% 
% controlMatrix = repmat([1e-10, 0, 0], [nSegment 1]); % play with numbers if you want
controlArray = reshape(controlMatrix', [4*nSegment, 1]);

% Converted to column vector, transpose needed to make the dimensions right
% InitialGuessPhase.state = stateArray;
% InitialGuessPhase.control = controlArray;
% InitialGuessPhase.timeSegment = time;
% InitialGuessPhase.nSegment = nSegment; % Remember that segment requires two nodes

InitialGuessPhase.state = stateMatrix;
InitialGuessPhase.control = controlMatrix;
InitialGuessPhase.timeSegment = time;
InitialGuessPhase.timeVariable = timeVariable;
InitialGuessPhase.timeDefect = timeDefect;
InitialGuessPhase.nSegment = nSegment; % Remember that segment requires two nodes
end