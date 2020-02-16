function StateMerged = mergePhase(State)

StateMerged = cell(1,1);

StateMerged{1}.state = [];
StateMerged{1}.control = [];
StateMerged{1}.timeSegment = [];
StateMerged{1}.nSegment = 0;
StateMerged{1}.initialConstraint = [];
StateMerged{1}.finalConstraint = [];
StateMerged{1}.slack = [];

nPhase = length(State);

for iPhase = 1:nPhase
	
	if iPhase == 1
		StateMerged{1}.state = [StateMerged{1}.state; State{iPhase}.state];
		StateMerged{1}.timeSegment = [StateMerged{1}.timeSegment; ...
			State{iPhase}.timeSegment];
		StateMerged{1}.initialConstraint = State{iPhase}.initialConstraint;
	else
		StateMerged{1}.state = [StateMerged{1}.state; State{iPhase}.state(8:end)];
		StateMerged{1}.timeSegment = [StateMerged{1}.timeSegment; ...
			State{iPhase}.timeSegment(2:end)];
	end
	
	if iPhase == nPhase
		StateMerged{1}.finalConstraint = State{iPhase}.finalConstraint;
	end
	
	StateMerged{1}.control = [StateMerged{1}.control; State{iPhase}.control];
	StateMerged{1}.nSegment = StateMerged{1}.nSegment + State{iPhase}.nSegment;

end