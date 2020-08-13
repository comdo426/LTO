function [indexSegmentOverlap, indexSegment] = getIndexSegment(StatePhase)
%GETINDEXSEGMENT - gets the index of the segment
%
%  Syntax:
%     [indexSegmentOverlap, indexSegment] = GETINDEXSEGMENT(StatePhase)
%
%  Description:
%     gets the indices of the segment.
%
%  Inputs:
%		StatePhase - structure of the state at a given phase
%
%  Outputs:
%     indexSegmentOverlap - index of the segment, both at the beginning of a
%     segment and at the end of a segment. Note that it only corresponds to the LGL
%     node placement scheme
%		indexSegment - index of the segment only at the beginning of a segment
%
%  See also: FMINCONCONSTRAINT
%
%   Author: Beom Park
%   Date: 18-Feb-2020; Last revision: 18-Feb-2020

indexSegmentOverlap = [];
indexSegment = [];

nSegment = StatePhase.nSegment;

for i = 1:nSegment
	indexSegmentOverlap = [indexSegmentOverlap, 4*(i-1)+1, 4*(i-1)+4]; %#ok<AGROW>
	indexSegment = [indexSegment, 4*(i-1)+1]; %#ok<AGROW>
end


end