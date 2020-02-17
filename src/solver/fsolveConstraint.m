function [F, dF] = fsolveConstraint(x, System, State, ...
	Spacecraft, Option, Collocation)
%FSOLVECONSTRAINT - computes the F and dF matrix for fsolve/newtonRaphson
%
%  Syntax:
%     [F, dF] = FSOLVECONSTRAINT(x, System, State, ...
%		Spacecraft, Option, Collocation)
%
%  Description:
%     Computes F and dF matrix for fsolve/newtonRaphson. 
%
%  Inputs:
%		x - parameter composed of state, control and slack variables
%
%  Outputs:
%     F - constraint column vector
%		dF - dF/dx matrix
%
%  See also: FMINCONCONSTRAINT
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 16-Feb-2020

nPhase = length(System);

%% number of continuity constraint

mContinuityTTL = 0;
isIniCon = ~isempty(State{1}.initialConstraint);
if isIniCon
	stateIniCon = State{1}.initialConstraint;
	nIniCon = length(stateIniCon);
	mContinuityTTL = mContinuityTTL + nIniCon;
end
isFinCon = ~isempty(State{end}.finalConstraint);
if isFinCon
	stateFinCon = State{end}.finalConstraint;
	nFinCon = length(stateFinCon);
	mContinuityTTL = mContinuityTTL + nFinCon;
end
mContinuityTTL = mContinuityTTL + 7*(nPhase-1);

phaseIni = zeros(7, nPhase);
phaseFin = zeros(7, nPhase);

%% number of parameter, number of constraints for each phase

[nx, ~] = getnx(State);
[mc, ~] = getmc(System, State, Spacecraft, Option);

%% dimension of the F, dF matrix

m = sum(mc) + mContinuityTTL; % row number of constraint
n = length(x); % column number of variables

F = nan(m, 1);
dF = zeros(m, n);

tstar = nan(nPhase,1);

for iPhase = 1:nPhase
	
	%% Information from the input
	
	isCR3BP = strcmp(System{iPhase}.dynamics(1), 'CR3BP');
	is2BP = strcmp(System{iPhase}.dynamics(1), '2BP');
	mu = System{iPhase}.parameter.mu;
	lstar = System{iPhase}.parameter.lstar;
	tstar(iPhase) = System{iPhase}.parameter.tstar;
	
	[t, nSegment, nNode, nState, nControl, mDefect, ~] = ...
		getPhaseStateInfo(State{iPhase});
	t = t + State{iPhase}.t0; % Add the offset from phasing(if needed)
	[thrustMaxND, ispND, g0ND] = getSpacecraftInfo(Spacecraft{iPhase});
	% nb: number of parameters before the i-th phase
	if iPhase == 1
		nb = 0;
		mb = 0;
	else
		nb = sum(nx(1:iPhase-1));
		mb = sum(mc(1:iPhase-1));
	end
	
	% index initialization
	cnt = 0;
	
	for i = 1:nSegment
		
		%% Defect constraint
		% define state at each node
		x1(:, 1) = x(28*(i-1)+1-7*cnt+nb:28*(i-1)+7-7*cnt+nb);
		x3(:, 1) = x(28*(i-1)+8-7*cnt+nb:28*(i-1)+14-7*cnt+nb);
		x5(:, 1) = x(28*(i-1)+15-7*cnt+nb:28*(i-1)+21-7*cnt+nb);
		x7(:, 1) = x(28*(i-1)+22-7*cnt+nb:28*(i-1)+28-7*cnt+nb);
		Y = [x1, x3, x5, x7];
		% define control at each segment
		u(:, 1) = x(nState + 3*(i-1)+1 +nb: nState + 3*i + nb);
		% time step at each segment
		dt = t(i+1) - t(i);
		% Defect constraint at each node
		
		
		if isCR3BP
			FDefect = defectConstraint(mu, Y, u, t(i), dt, Collocation, ispND, g0ND);
			[dFx, dFu] = dFDefect(mu, Y, u, t(i), dt, Collocation, ispND, g0ND);
		end
		
		if is2BP
			FDefect = defectConstraint2BP(mu, Y, u, t(i), dt, Collocation, ispND, g0ND);
			[dFx, dFu] = dFDefect2BP(mu, Y, u, t(i), dt, Collocation, ispND, g0ND);
		end
		
		F(21*(i-1)+1 + mb: 21*i + mb, 1) = FDefect;
		dF(21*(i-1)+1+mb: 21*i+mb, 28*(i-1)+1-7*cnt+nb: 28*(i-1)+28-7*cnt+nb) = dFx;
		dF(21*(i-1)+1+mb: 21*i+mb, nState+3*(i-1)+1+nb: nState+3*i+nb) = dFu;
		
		% slack variable dF at each node
		
		%% Thrust slack constraint
		% slack variable at each segment
		lambda = x(nState + nControl + i + nb);
		mThrust = nSegment;
		nThrust = nSegment;
		% slack variable constraint at each node
		% 		FPhase(nDefect + i, 1) = u(1) - thrustMaxND*sin(lambda)^2;
		% 		dFPhase(nDefect+i, nState+3*(i-1)+1) = 1;
		% 		dFPhase(nDefect+i, nState+nControl+i) = -2*thrustMaxND*sin(lambda)*cos(lambda);
		F(mDefect + i + mb, 1) = u(1) - thrustMaxND*sin(lambda)^2;
		dF(mDefect + i + mb, nState + 3*(i-1)+1 + nb) = 1;
		dF(mDefect + i + mb, nState + nControl + i +nb) = ...
			-2*thrustMaxND*sin(lambda)*cos(lambda);
		
		if i == 1
			phaseIni(:, iPhase) = x1;
		end
		if i == nSegment
			phaseFin(:, iPhase) = x7;
		end
		% increment on index
		cnt = cnt + 1;
	end % iSegment for loop
	
	
	%% Additional constraints
	if ~isempty(Option.AddCon{iPhase, 1})
		
		for iNode = 1:nNode
			
			sigma1 = x(nState + nControl + nThrust + iNode + nb);
			nSigma1 = nNode;
			mSigma1 = nNode;
			sigma2 = x(nState + nControl + nThrust + nSigma1 + iNode +nb);
			
			minAlt1 = Option.AddCon{iPhase,1}{5}/lstar;
			minAlt2 = Option.AddCon{iPhase,2}{5}/lstar;
			pos = x(7*(iNode-1)+1+nb:7*(iNode-1)+3+nb);
			alt1 = norm(pos' - [-mu, 0, 0]);
			alt2 = norm(pos' - [1-mu, 0, 0]);
			
			F(mDefect + mThrust + iNode + mb, 1) = ...
				alt1^2*sin(sigma1) - minAlt1^2;
			F(mDefect + mThrust + mSigma1 + iNode + mb, 1) = ...
				alt2^2*sin(sigma2) - minAlt2^2;
			
			dF(mDefect + mThrust + iNode + mb, 7*(iNode-1)+1 + nb:7*(iNode-1)+3 + nb) = ...
				2*(pos' - [-mu, 0, 0])*sin(sigma1);
			dF(mDefect + mThrust + iNode + mb, nState+nControl+nThrust+iNode +nb) = ...
				alt1^2*cos(sigma1);
			dF(mDefect + mThrust + mSigma1 + iNode + mb, 7*(iNode-1)+1 +nb:7*(iNode-1)+3 + nb) = ...
				2*(pos' - [1-mu, 0, 0])*sin(sigma2);
			dF(mDefect + mThrust + mSigma1 + iNode + mb, nState+nControl+nThrust+nSigma1+iNode +nb) = ...
				alt2^2*cos(sigma2);
			
		end
	end
	
	%% Continuity constraint
	
	mt = sum(mc); % number of constraints over the phases
	
	% F
	if iPhase == 1
		F(mt+1: mt+nIniCon) = phaseIni(:, 1) - stateIniCon;
	else
		isSame = isequaln(System{iPhase}.dynamics, System{iPhase-1}.dynamics);
		if isSame % Share same dynamics
			F(mt+nIniCon+7*(iPhase-2)+1: mt+nIniCon+7*(iPhase-1)) = ...
				phaseIni(:,iPhase) - phaseFin(:,iPhase-1);
		else % different dynamics: need to rotate
			JDAdd = State{iPhase-1}.timeSegment(end)*tstar(iPhase-1)/60/60/24;
			phaseFinNew = interPhaseCont(phaseFin(:,iPhase-1), System{iPhase-1}, ...
				System{iPhase}, State{iPhase-1}.JD0+JDAdd);
			F(mt+nIniCon+7*(iPhase-2)+1: mt+nIniCon+7*(iPhase-1)) = ...
				phaseIni(:,iPhase) - phaseFinNew;
		end
	end
	if iPhase == nPhase
		F(mt+nIniCon+7*(iPhase-1)+1: mt+nIniCon+7*(iPhase-1)+nFinCon) = ...
			phaseFin(1:nFinCon,iPhase) - stateFinCon;
	end
	
	% dF
	if iPhase == 1
		dF(mt+1: mt+nIniCon, 1:nIniCon) = eye(nIniCon);
		if nPhase > 1
		isSame = isequaln(System{iPhase+1}.dynamics, System{iPhase}.dynamics);
			if isSame
				dF(mt+nIniCon+1: mt+nIniCon+7, (nState-7)+1:nState) = -eye(7);
			else
				JDAdd = State{iPhase}.timeSegment(end)*tstar(iPhase)/60/60/24;
				JD = State{iPhase}.JD0 + JDAdd;
				dF(mt+nIniCon+1: mt+nIniCon+7, (nState-7)+1:nState) = dFInterPhaseCont(...
					phaseFin(:,iPhase), System{iPhase}, System{iPhase+1}, JD);
			end
		end
	else
		dF(mt+nIniCon+7*(iPhase-2)+1: mt+nIniCon+7*(iPhase-1), ...
			1+nb: 7+nb) = eye(7);
	end
	
	if iPhase == nPhase
		dF(mt+nIniCon+7*(iPhase-1)+1: mt+nIniCon+7*(iPhase-1)+nFinCon, ...
			(nState-7)+1+nb: (nState-7)+nFinCon+nb) = eye(nFinCon);
	else
		if iPhase > 1
			dF(mt+nIniCon+7*(iPhase-1)+1: mt+nIniCon+7*iPhase, ...
				(nState-7)+1+nb: nState+nb) = -eye(7);
		end
	end
	
end % iPhase for loop

save('TEST')

end