% Caution, only works for the altitude constraint

function xShrink = deleteSlackVariable(x, State, Option, feasible, optimize)

nPhase = length(State);

xShrink = [];

if feasible
for iPhase = 1:nPhase
	if iPhase == 1
		nx(iPhase) = 0;
	else
		nSegment = State{iPhase-1}.nSegment; % number before the phase
		if ~isempty(Option.AddCon{iPhase-1, 1})
			nx(iPhase) = nx(iPhase-1) + 21*nSegment + 7 + 3*nSegment + nSegment + (3*nSegment + 1)*2;
		else
			nx(iPhase) = nx(iPhase-1) + 21*nSegment + 7 + 3*nSegment + nSegment; % with only thrust
		end
	end
end

for iPhase = 1:nPhase
	nSegment = State{iPhase}.nSegment;
	nSlackThrust = nSegment;
	if ~isempty(Option.AddCon{iPhase, 1})
		nSlackAltitude = (3*nSegment + 1)*2;
	else
		nSlackAltitude = 0;
	end
	nNode = 3*nSegment + 1;
	nState = 7*nNode;
	nControl = 3*nSegment;
	nState+nControl+1+nx(iPhase);
	nState+nControl+nSlackThrust+nSlackAltitude+nx(iPhase);
	x(nState+nControl+1+nx(iPhase):nState+nControl+nSlackThrust+nSlackAltitude+nx(iPhase)) = [];
	nx = nx-nSlackThrust-nSlackAltitude;
% 	xShrink = [xShrink; x(nx(iPhase)+1:nx(iPhase)+nState+nControl)];
end
end

if optimize
for iPhase = 1:nPhase
	if iPhase == 1
		nx(iPhase) = 0;
	else
		nSegment = State{iPhase-1}.nSegment; % number before the phase
		if ~isempty(Option.AddCon{iPhase, 1})
			nx(iPhase) = 21*nSegment + 7 + 3*nSegment + (3*nSegment + 1)*2;
		else
			nx(iPhase) = 21*nSegment + 7 + 3*nSegment; % with only thrust
		end
	end
end

for iPhase = 1:nPhase
	nSegment = State{iPhase}.nSegment;
	nSlackThrust = 0;
	if ~isempty(Option.AddCon{iPhase, 1})
		nSlackAltitude = (3*nSegment + 1)*2;
	else
		nSlackAltitude = 0;
	end
	nNode = 3*nSegment + 1;
	nState = 7*nNode;
	nControl = 3*nSegment;
	x(nState+nControl+1+nx(iPhase):nState+nControl+nSlackThrust+nSlackAltitude+nx(iPhase)) = [];
	nx = nx-nSlackThrust-nSlackAltitude;
end
	

end

xShrink = x;
end