function State = updateState(State, x)

nPhase = length(State);

for iPhase = 1:nPhase
	
	State{iPhase}.slack = [];
end

nx = getnx(State);

for iPhase = 1:nPhase
	[~, ~, ~, nState, nControl, ~, ~] = getPhaseStateInfo(State{iPhase});
	
	if iPhase == 1
		nb = 0;
	else
		nb = sum(nx(1:iPhase-1));
	end
	State{iPhase}.state = x(1+nb:nState+nb);
	State{iPhase}.control = x(1+nState+nb:nState+nControl+nb);
	
	
end





end