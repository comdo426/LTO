function State = LTOMain(System, State, Spacecraft, Option)

%% set boolean variables

isFeasible = LToptget(Option.LTO, 'FeasibilitySolver');
isOptimize = LToptget(Option.LTO, 'Optimizer');
isMesh = LToptget(Option.LTO, 'MeshRefinement');

%% Collocation matrices

[phi, phiPrime, phiMeshAdd] = LGL_7th_coefficient;
Collocation.phi = phi;
Collocation.phiPrime = phiPrime;
Collocation.phiMeshAdd = phiMeshAdd;
Collocation.tau = [-1, -0.830223896278567, -0.468848793470714, 0, ...
	0.468848793470714, 0.830223896278567, 1];
tau = Collocation.tau;
Collocation.tauRatio = (tau - ones(1, 7)*tau(1))/2;

%% feasible

iWhile = 1;

while isFeasible && ~Option.doneFeasible
	fprintf('Feasible: step no. %d\n', iWhile);
	Problem = setProblemfsolve(System, State, Spacecraft, Option, ...
		Collocation);
% 	save('Problem')
% 	break
	if ~isempty(Option.newton) % newton-raphson method
		feasibleWithSlack = newtonRaphson(Problem);
		feasibleVec = deleteSlackVariable(feasibleWithSlack, State, Option, 1, 0);
		State = updateState(State, feasibleVec);
      closestEncounter(System, State);
		save('TEST2')
	else % try fsolve
		feasibleWithSlack = fsolve(Problem);
		feasibleVec = deleteSlackVariable(feasibleWithSlack, State);
		State = updateState(State, feasibleVec);
      closestEncounter(System, State);
	end % newton-raphson if loop
	
	if Option.plot.feasible{1}
		[earthPlot, moonPlot, lpPlot] = drawThrustArc(Option.plot.feasible{2}(1,iWhile), Option.plot.feasible{3}, ...
			System, State, Spacecraft);
      [initialPlot, interPhase, finalPlot] = drawEndPoints(Option.plot.feasible{2}(1,iWhile), System, State);
		drawThrustHistory(Option.plot.feasible{2}(2,iWhile), ...
			System, State, Spacecraft);
% 		figure(Option.plot.feasible{2}(1,iWhile))
% 		legend([initialPlot, interPhase, finalPlot, ...
% 			 moonPlot, lpPlot], {'Initial', 'InterPhase', 'Final', ...
% 		 'Moon', 'L_i'}, 'fontsize', 12);
		figure(Option.plot.feasible{2}(1,iWhile)+100)
		legend([initialPlot, interPhase, finalPlot], {'Initial', 'InterPhase', ...
			'Final'}, 'fontsize', 13);
	end
	
	
	if isMesh
		[State, Option] = ...
			CEPMeshRefineMultiPhase(System, State, Spacecraft, Option);
		if Option.doneMesh
			Option.doneFeasible = 1;
		end
	else
		Option.doneFeasible = 1;
	end % isMesh if loop
	iWhile = iWhile + 1;
end % isFeasible, ~Option.doneFeasible while loop

prompt = 'Feasible done! Wish to continue?';
C = input(prompt);


%% Optimize

jWhile = 1;

while isOptimize && ~Option.doneOptimize
	fprintf('Optimize: step no. %d\n', jWhile);
	Problem = setProblemfmincon(System, State, Spacecraft, Option);
	optimizedWithSlack = fmincon(Problem);
	optimizedVec = deleteSlackVariable(optimizedWithSlack, State, Option, 0, 1);
	State = updateState(State, optimizedVec);
      closestEncounter(System, State);
	
	if Option.plot.optimize{1}
		[~, moonPlot, lpPlot] = drawThrustArc(Option.plot.optimize{2}(1,jWhile), Option.plot.optimize{3}, ...
			System, State, Spacecraft);
      [initialPlot, interPhase, finalPlot] = drawEndPoints(Option.plot.optimize{2}(1,jWhile), System, State);
		drawThrustHistory(Option.plot.optimize{2}(2,jWhile), ...
			System, State, Spacecraft);
% 		figure(Option.plot.optimize{2}(1,jWhile))
% 		legend([initialPlot, interPhase, finalPlot, ...
% 			moonPlot, lpPlot], {'Initial', 'InterPhase', 'Final', ...
% 			'Earth', 'Moon', 'L_i'}, 'fontsize', 12);
		figure(Option.plot.optimize{2}(1,jWhile)+100)
		legend([initialPlot, interPhase, finalPlot], {'Initial', 'InterPhase', ...
			'Final'}, 'fontsize', 13);
	end
	
	if isMesh
		[State, Option] = ...
			CEPMeshRefineMultiPhase(System, State, Spacecraft, Option);
		if Option.doneMesh
			Option.doneOptimize = 1;
		end
	else
		Option.doneOptimize = 1;
	end
	jWhile = jWhile + 1;
end % isOptimize, ~Option.doneOptimize while loop

end