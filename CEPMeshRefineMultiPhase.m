function [State, Option] = ...
	CEPMeshRefineMultiPhase(System, State, Spacecraft, Option)

nPhase = length(State);
doneMesh = zeros(nPhase,1);

for j = 1:nPhase
	
	j

	if j < 2 % Ad hoc version
	[State{j}, doneMesh(j,1)] = ...
		CEPMeshRefineSinglePhase(System{j}, State{j}, Spacecraft{j}, Option);
	doneMesh(2) = 1;
	end
	
end

doneMesh

if all(doneMesh == 1)
	Option.doneMesh = 1;
else
	Option.doneMesh = 0;
end
   
end