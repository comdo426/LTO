function n = getAltitudeConstraintNo(~, StatePhase, ...
	~, ~)

[~, ~, nNode, ~, ~, ~, ~] = getPhaseStateInfo(StatePhase);

n = nNode;


end