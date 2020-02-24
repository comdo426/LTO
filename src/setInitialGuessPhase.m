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
stateMatrix = nan(3*nSegment+1, 7);
for iRev = 1:intRev
	if iRev < intRev
		stateMatrix(3*s*(iRev-1)+1:3*s*iRev, :) = stateOneRev;
	elseif isInt
		stateMatrix(3*s*(iRev-1)+1:3*s*iRev+1, :) = ...
			[stateOneRev; [stateInitial, 1]];
	else
		stateMatrix(3*s*(iRev-1)+1:3*s*iRev, :) = stateOneRev;
		stateMatrix(3*s*iRev+1:3*nSegment+1, :) = stateExtra;
	end % if loop
end % iRev for loop
stateArray = reshape(stateMatrix', [7*(3*nSegment + 1), 1]);

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

% control formulation
controlMatrix = repmat([1e-10, 0, 0], [nSegment 1]); % play with numbers if you want
controlArray = reshape(controlMatrix', [3*nSegment, 1]);

% Converted to column vector, transpose needed to make the dimensions right
InitialGuessPhase.state = stateArray;
InitialGuessPhase.control = controlArray;
InitialGuessPhase.timeSegment = time;
InitialGuessPhase.nSegment = nSegment; % Remember that segment requires two nodes

end