function [f g] = fminconObjective(x, iFinalMass)
f = -x(iFinalMass);
if nargout > 1
	g = zeros(length(x), 1);
	g(iFinalMass) = -1;
end
end