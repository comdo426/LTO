function m = getAltitudeConstraintNo(~, StatePhase, ~, ~)
%GETALTITUDECONSTRAINTNO - gets the number of constraints for the altitude cont.
%
%  Syntax:
%     M = getAltitudeConstraintNo(~, StatePhase, ~, ~)
%	
%  Description:
%     gets the number of constraints for the each iPhase. 
%
%  Outputs:
%     m - scalar, number of altitude constraints for iPhase
%
%  See also: SETPROBLEMFSOLVE, SETPROBLEMFMINCON, GETNX
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 16-Feb-2020

[~, ~, nNode, ~, ~, ~, ~] = getPhaseStateInfo(StatePhase);

m = nNode;

end