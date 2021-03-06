function Feasible = newtonRaphson(Problem)
%NEWTONRAPHSON - NewtonRaphson method to solve for zero
%
%  Syntax:
%     Feasible = NEWTONRAPHSON(Problem)
%
%  Description:
%     Uses NewtonRaphson method to solve for zero of the objective function
%
%  Inputs:
%		Problem - structure that contains the following information
%			x0 - initial guess
%			options - structure with
%				fcnTolerance - the tolerance to stop the iteration
%				maxIteration - maximum number of iteration
%			objective - function handle, it is set to be fsolveConstraint
%
%  Outputs:
%     Feasible - column vector of converged state, control and also slack variables.
%     Note that you need another function(deleteSlackVariable) to run it again
%     after mesh refinement
%
%  See also: SETPROBLEMFSOLVE, FSOLVECONSTRAINT, DELETESLACKVARIABLE
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 16-Feb-2020


tol = Problem.options.fcnTolerance;
iNewtonMax = Problem.options.maxIteration;

iNewton = 1;
normF = 1000;
dx = zeros(length(Problem.x0), 1);
while normF > tol
	if iNewton == 1
		x = Problem.x0;
	end
	if iNewton > iNewtonMax
		break;
	end
	[F, dF] = Problem.objective(x);
	normF = norm(F);
	fprintf('NewtonRaphson step no. %d: error norm is %e\n', iNewton, normF);
	
	dx = - dF'*((dF*dF')\F);
	x = x+dx;
	
	iNewton = iNewton + 1;
end
if iNewton
x = x - dx; % return to the original solution, 

Feasible = x;

end