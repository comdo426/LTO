function [t, nSegment, nNode, nState, nControl, nDefect, nSlack] = ...
	getPhaseStateInfo(StatePhase)

t = StatePhase.timeSegment;
nSegment = StatePhase.nSegment;
nNode = 3*nSegment + 1;
nState = 7*nNode;
nControl = 3*nSegment;
nDefect = 21*nSegment;
nSlack = length(StatePhase.slack);

end