function [stateMat, stateSegmentMat, controlMat] = ...
	getStateControlMat(StatePhase)

state = StatePhase.state;
control = StatePhase.control;
[~, nSegment, nNode, nState, nControl, ~, ~] = getPhaseStateInfo(StatePhase);

stateMat = reshape(state, [7, nNode])'; %Transpose

stateSegmentMat = zeros(nSegment+1,7);
stateSegmentMat(1,:) = stateMat(1,:);
for i = 1:nSegment
	stateSegmentMat(i+1,:) = stateMat(3*i+1, :);
end

controlMat = reshape(control, [3, nSegment])'; %Transpose
	
end