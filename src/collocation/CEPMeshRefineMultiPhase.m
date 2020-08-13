function [State, Option] = ...
	CEPMeshRefineMultiPhase(System, State, Spacecraft, Option)
%CEPMESHREFINEMULTIPHASE - refines mesh with CEP(Control Explicit Propagation)
%
%  Syntax:
%     [State, Option] = ...
%		CEPMESHREFINEMULTIPHASE(System, State, Spacecraft, Option)
%
%  Description:
%     refines mesh with CEP method.
%
%  Inputs:
%		Option - contains the information about the mesh refinement settings
%		Example:
%			Option.doneFeasible = 0;
%			Option.doneOptimize = 0;
%			Option.doneMesh = 0;
%			Option.removeMesh = 0;
%			Option.meshTolerance = 1e-12;
%
%  Outputs:
%     Option -
%			doneMesh - boolean that denotes the end of mesh refine
%
%	See also: CEPMESHREFINESINGLEPHASE
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 17-Feb-2020

nPhase = length(State);
doneMesh = zeros(nPhase,1);

for j = 1:nPhase
	
	cprintf(-[0 1 1], 'Phase no. %d, mesh refinement\n', j);
	[State{j}, doneMesh(j,1)] = ...
		CEPMeshRefineSinglePhase(System{j}, State{j}, Spacecraft{j}, Option);
	
end

if Option.removeMesh
	if all(doneMesh == 1)
		Option.removeMesh = false;
	else
		Option.removeMesh = true;
	end
else
	if all(doneMesh == 1)
		Option.doneMesh = true;
	else
		Option.doneMesh = false;
	end
end

end