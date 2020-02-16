function Problem = setProblemfsolve(System, State, Spacecraft, Option, ...
	Collocation)

nPhase = length(State);

%% initialize

Problem.x0 = [];

for iPhase = 1:nPhase
	[~, nSegment, nNode, ~, ~, ~, ~] = getPhaseStateInfo(State{iPhase});
	Tmax = Spacecraft{iPhase,1}.thrustMaxND;
	lambda = nan(nSegment, 1);
	[stateMat, ~, controlMat] = getStateControlMat(State{iPhase});
	for iSegment = 1:nSegment
		T = controlMat(iSegment, 1);
		lambda(iSegment, 1) = asin(sqrt(T/Tmax));
	end
	Problem.x0 = [Problem.x0; State{iPhase}.state; ...
		State{iPhase}.control; lambda];
	State{iPhase}.slack = lambda;
	
	if ~isempty(Option.AddCon{iPhase, 1})
		mnAddCon = size(Option.AddCon);
		nAddCon = mnAddCon(1, 2);
		
		for iAddCon = 1:nAddCon
			Con = Option.AddCon{iPhase, iAddCon};
			isSlack = strcmp(Con{2}, 'ineq'); % If slack, we need to augment the x0
			if isSlack
				c = str2func(strcat('get',Con{3},'Constraint'));
				sigma = c(System{iPhase}, State{iPhase}, Spacecraft{iPhase}, Con);
				Problem.x0 = [Problem.x0; sigma];
				State{iPhase}.slack = [State{iPhase}.slack; sigma];
			end
		end
	end
	
end


%% set options
if ~isempty(Option.newton)
	Problem.options = Option.newton;
else
	Problem.options = Option.fsolve;
end
Problem.solver = 'fsolve';

save('State', 'State');
save('Problem', 'Problem');

%% set functions

[phi, phiPrime, phiMeshAdd] = LGL_7th_coefficient;
Collocation.phi = phi;
Collocation.phiPrime = phiPrime;
Collocation.phiMeshAdd = phiMeshAdd;
Collocation.tau = [-1, -0.830223896278567, -0.468848793470714, 0, ...
	0.468848793470714, 0.830223896278567, 1];
tau = Collocation.tau;
Collocation.tauRatio = (tau - ones(1, 7)*tau(1))/2;

Problem.objective = @(x) fsolveConstraint(x, System, State, Spacecraft, ...
	Option, Collocation); % working on it!

end