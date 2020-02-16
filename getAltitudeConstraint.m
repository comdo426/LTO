function sigma = getAltitudeConstraint(SystemPhase, StatePhase, ...
	~, Con)
% Altitude constraint. Note that constraint is set to be x^2 + y^2 + z^2 >= h^2

mu = SystemPhase.parameter.mu;
lstar = SystemPhase.parameter.lstar;

[~, ~, nNode, ~, ~, ~, ~] = getPhaseStateInfo(StatePhase);
[stateMat, ~, ~] = getStateControlMat(StatePhase);

primaryNo = Con{4};
minAlt = Con{5}/lstar;

isCR3BP = strcmp(SystemPhase.dynamics{1}, 'CR3BP');

sigma = nan(nNode, 1);

for iNode = 1:nNode
	if isCR3BP
		if primaryNo == 1
			alt = norm(stateMat(iNode, 1:3) - [-mu, 0, 0]);
		else
			if primaryNo == 2
				alt = norm(stateMat(iNode, 1:3) - [1-mu, 0, 0]);
			else
				fprintf('primaryNo input Wrong');
			end
		end
	else
		fprintf('Dynamics other than CR3BP not supported');
	end
	
	if alt >= minAlt
		sigma(iNode, 1) = asin(minAlt^2/alt^2);
	else
		fprintf('Altitude constraint violation: %0.5e, choose orbit with higher altitude\n', ...
			alt);
	end
end

end