% function Feasible = newtonRaphson(System, InitialGuess, Spacecraft, ...
% 	Option, Collocation, Problem)
function Feasible = newtonRaphson(Problem)

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