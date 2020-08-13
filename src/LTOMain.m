function State = LTOMain(System, State, Spacecraft, Option)
%LTOMAIN - solves the low-thrust optimization problem
%
%  Syntax:
%     State = LTOMAIN(System, State, Spacecraft, Option)
%
%  Description:
%     Solves the LTO problem with LGL 7th order collocation scheme
%
%  Inputs:
%		System: cell for the dynamics and parameter of each phase
%		State: cell for the state, control, constraints of each phase
%		Spacecraft: cell for the dimensional/non-dimensional spacecraft specs
%		Option: cell for the options of integration, and solvers
%
%  Outputs:
%     State - cell for the state, control of the converged solution
%
%   Author: Beom Park
%   Date: 16-Feb-2020; Last revision: 01-Mar-2020

%% Settings before going into the algorithm

% boolean variables
isFeasible = LToptget(Option.LTO, 'FeasibilitySolver');
isOptimize = LToptget(Option.LTO, 'Optimizer');
isMesh = LToptget(Option.LTO, 'MeshRefinement');

% Collocation structure
Collocation = setCollocation;

%% feasible

iWhile = 1;

while isFeasible && ~Option.doneFeasible
	cprintf(-[1, 0, 0], 'Feasible: step no. %d\n', iWhile);
	[Problem, State] = setProblemfsolve(System, State, Spacecraft, Option, ...
		Collocation);
	if ~isempty(Option.newton) % newton-raphson method
		feasibleWithSlack = newtonRaphson(Problem);
	else % try fsolve
		feasibleWithSlack = fsolve(Problem);
	end % newton-raphson if loop
	feasibleVec = deleteSlackVariable(feasibleWithSlack, State, Option);
	State = updateState(State, feasibleVec);
	closestEncounter(System, State);
	
	if Option.plot.feasible{1}
		[~, ~, ~] = drawThrustArc(Option.plot.feasible{2}(1,iWhile), Option.plot.feasible{3}, ...
			System, State, Spacecraft);
		[initialPlot, interPhase, finalPlot] = drawEndPoints(Option.plot.feasible{2}(1,iWhile), System, State);
		drawThrustHistory(Option.plot.feasible{2}(2,iWhile), ...
			System, State, Spacecraft);
		% plot Earth, Moon, etc
		% 		figure(Option.plot.feasible{2}(1,iWhile))
		% 		legend([initialPlot, interPhase, finalPlot, ...
		% 			 moonPlot, lpPlot], {'Initial', 'InterPhase', 'Final', ...
		% 		 'Moon', 'L_i'}, 'fontsize', 12);
		figure(Option.plot.feasible{2}(1,iWhile)+100)
		legend([initialPlot, interPhase, finalPlot], {'Initial', 'InterPhase', ...
			'Final'}, 'fontsize', 13);
	end % plot for feasible if loop
		
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

cprintf(-[1, 0, 0], 'Feasible done!\n');

%% Optimize

jWhile = 1;

while isOptimize && ~Option.doneOptimize
	cprintf(-[1, 0, 0], 'Optimize: step no. %d\n', jWhile);
	[x0, funcs, options, State] = setProblemipopt(System, State, Spacecraft, Option, Collocation);
	[optimizedWithSlack, ~] = ipopt_auxdata(x0, funcs, options);
	optimizedVec = deleteSlackVariable(optimizedWithSlack, State, Option);
	State = updateState(State, optimizedVec);
	closestEncounter(System, State);
	
	if Option.plot.optimize{1}
		drawThrustArc(Option.plot.optimize{2}(1,jWhile), Option.plot.optimize{3}, ...
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