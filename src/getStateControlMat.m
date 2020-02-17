function [stateMat, stateSegmentMat, controlMat] = ...
	getStateControlMat(StatePhase)
%GETSTATECONTROLMAT - gets the state/control matrices of State{iPhase}
%
%  Syntax:
%     [stateMat, stateSegmentMat, controlMat] = ...
%		GETSTATECONTROLMAT(StatePhase)
%
%  Description:
%     converts the given State{iPhase} into a matrix form. This is helpful when
%     we have to draw plots.
%
%  Inputs:
%		StatePhase, or State{iPhase}
%
%  Outputs:
%     stateMat - state matrix for all the nodes in it
%		stateSegmentMat - state matrix for only the segments
%		controlMat - control matrix for the segments
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 16-Feb-2020

state = StatePhase.state;
control = StatePhase.control;
[~, nSegment, nNode, ~, ~, ~, ~] = getPhaseStateInfo(StatePhase);

stateMat = reshape(state, [7, nNode])'; %Transpose to match the dimensions

stateSegmentMat = zeros(nSegment+1,7);
stateSegmentMat(1,:) = stateMat(1,:);
for i = 1:nSegment
	stateSegmentMat(i+1,:) = stateMat(3*i+1, :);
end

controlMat = reshape(control, [3, nSegment])'; %Transpose
	
end