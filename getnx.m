function [nx, nxOpt] = getnx(State)


nPhase = length(State);
nx = zeros(nPhase, 1);
nxOpt = zeros(nPhase, 1);

for iPhase = 1:nPhase
	[~, ~, ~, nState, nControl, ~, nSlack] = getPhaseStateInfo(State{iPhase});
	nx(iPhase) = nState + nControl + nSlack;
	nxOpt(iPhase) = nState + nControl + nSlack;
	% works for altitude constraints
end

end