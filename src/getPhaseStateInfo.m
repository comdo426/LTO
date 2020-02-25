function [t, nSegment, nNode, nState, nControl, mDefect, nSlack] = ...
	getPhaseStateInfo(StatePhase)
%GETPHASESTATEINFO - gets the information per state in a certain phase
%
%  Syntax:
%     [t, nSegment, nNode, nState, nControl, mDefect, nSlack] = ...
%		GETPHASESTATEINFO(StatePhase)
%
%  Description:
%     gets the significant numbers of the State{iPhase} structure
%
%  Outputs:
%     t: time vector for the segment
%		nSegment: number of segment
%		nNode: number of nodes
%		nState: number of states
%		nControl: number of control variables
%		mDefect: number of defect constraints. m is used to denote it's the
%		constraint(row, from m-n notation)
%		nSlack: number of slack variables BESIDES the default thrust slack
%		variables
%
%  See also: SETINITIALGUESS(for the State{iPhase}.slack initialization)
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 16-Feb-2020

t = StatePhase.timeSegment;
nSegment = StatePhase.nSegment;
nNode = 4*nSegment;
nState = 7*nNode;
nControl = 3*nSegment;
mDefect = 21*nSegment;
nSlack = length(StatePhase.slack);

end