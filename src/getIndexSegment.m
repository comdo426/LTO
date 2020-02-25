function indexSegment = getIndexSegment(StatePhase)

indexSegment = [];

nSegment = StatePhase.nSegment;

for i = 1:nSegment
	indexSegment = [indexSegment, 4*(i-1)+1, 4*(i-1)+4]; %#ok<AGROW>
end


end