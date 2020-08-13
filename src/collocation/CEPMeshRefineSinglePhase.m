function [StatePhase, doneMeshPhase] = ...
	CEPMeshRefineSinglePhase(SystemPhase, StatePhase, SpacecraftPhase, Option)
%CEPMESHREFINESINGLEPHASE - refines mesh with CEP single phase
%
%  Syntax:
%     [StatePhase, doneMeshPhase] = ...
%		CEPMESHREFINESINGLEPHASE(SystemPhase, StatePhase, SpacecraftPhase, Option)
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
%     StatePhase - structure for the refined State for each phase
%		doneMeshPhase - boolean variable denoting the end of mesh refinement
%
%	See also: CEPMESHREFINEMULTIPHASE, CR3BPLT_CSI_INERTIAL, CONICLT
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 01-Mar-2020

% Define the propagator depending on the dynamics
isCR3BP = strcmp(SystemPhase.dynamics{1}, 'CR3BP');
is2BP = strcmp(SystemPhase.dynamics{1}, '2BP');
if isCR3BP
	prop = str2func('CR3BPLT');
end
if is2BP
	prop = str2func('conicLT');
end

% Set necessary numbers
t0 = StatePhase.t0;
opts = odeset('reltol', 3e-14, 'abstol', 1e-20);
tol = Option.meshTolerance;
mu = SystemPhase.parameter.mu;
[t, s, ~, nState, nControl, ~, ~] = getPhaseStateInfo(StatePhase);
[~, ispND, g0ND] = getSpacecraftInfo(SpacecraftPhase);

% Put state/control into a single vector x
% x = [StatePhase.state; StatePhase.control];
x = [reshape(StatePhase.state', [nState, 1]); reshape(StatePhase.control', ...
	[nControl, 1])];

% Put state/control into a matrix x_s, u
% [~, x_s, u] = getStateControlMat(StatePhase);
indexSegment = getIndexSegment(StatePhase);
x_s = StatePhase.state(indexSegment, :);
u = StatePhase.control;

% Set collocation variable/matrices
Collocation = setCollocation;
invA = Collocation.invA;
BAdd = Collocation.BAdd;

tauRatio = Collocation.tauRatio;
% phiMeshAdd = Collocation.phiMeshAdd;

if Option.removeMesh
	
	%% Remove
	
	k = 1;
	rN = []; % short word for remove number
	skipSwitch = 0;
	normER = nan(s-1, 1);
	for i = 1:s-1
		if skipSwitch ~= 1
			% First propagation
			[~, state] = ode113(@(t,y) prop(t,y,u(i,:),mu, ...
				ispND,g0ND,t0), [t(i), t(i+1)], x_s(i,:), opts);
			% Second propagation
			[~, state] = ode113(@(t,y) prop(t,y,u(i,:),mu, ...
				ispND,g0ND,t0), [t(i+1), t(i+2)], state(end,:), opts);
			xProp = state(end, :);
			E = xProp - x_s(i+2,:);
			normER(i) = norm(E);
			if normER(i) < tol/100
				rN(k) = i;
				fprintf('Remove: Segment no. %d, the error is %0.5e\n', ...
					rN(k), normER(i));
				k = k+1;
				skipSwitch = 1;
			end
		else
			skipSwitch = 0;
		end
	end
	
	if ~isempty(rN)
		nRemove = length(rN);
		xNew = [];
		uNew = [];
		for j = 1:nRemove
			% Find the x state at the beginning of the current segment and
			% end of the segment
			xj1 = x(28*(rN(j)-1)+1:28*(rN(j)-1)+7);
			xj2 = x(28*(rN(j)+1)+1:28*(rN(j)+1)+7);
			xR(:, j) = [xj1+(xj2-xj1)*tauRatio(3); xj1+(xj2-xj1)*tauRatio(5)];
			if nRemove == 1
				xNew = [x(1:21*(rN(j)-1)+7); xR; x(21*(rN(j)+1)+nState)];
				uNew = [x(nState+1: nState+4*(rN(j)-1)+4); x(nState+4*nR(j)+5:end)];
			else
				if j < nRemove
					if j == 1
						xNew = [x(1:28*(rN(j)-1)+7); xR];
						uNew = [x(nState+1: nState+4*(rN(j)-1)+4)];
					else
						xNew = [xNew; x(21*(nR(j-1)+1)+1: 21*(nR(j)-1)+7); xR];
						uNew = [uNew; x(nState+4*(nR(j-1)+1)+1: nState+4*(nR(j)-1)+4)];
					end
				else
					xNew = [xNew; x(28*(nR(j-1)+1)+1: 28*(nR(j)-1)+7); xR; ...
						x(28*(nR(j)+1)+1: nState)];
					uNew = [uNew; x(nState+4*(nR(j-1)+1)+1: nState+4*(nR(j)-1)+4); ...
						x(nState+4*nR(j)+5:end)];
				end
			end
		end
		StatePhase.state = xNew;
		StatePhase.control = uNew;
		StatePhase.nSegment = length(t)-1;
		doneMeshPhase = false;
	else
		doneMeshPhase = true;
	end
	
else
	
	%% Add
	k = 1;
	aN = [];
	normEA = nan(s, 1);
	for i = 1:s
		[~, state] = ode113(@(t,y) prop(t,y,u(i,:),mu,ispND, ...
			g0ND, t0), [t(i), t(i+1)], x_s(2*(i-1)+1,:), opts);
		xProp = state(end,:);
		E = xProp - x_s(2*i, :);
		normEA(i) = norm(E);
		if normEA(i) > tol
			aN(k) = i;
			fprintf('Add: Segment no. %d, the error is %0.5e\n', ...
				aN(k), normEA(i));
			k = k+1;
		end
	end
	
	if ~isempty(aN)
		cprintf(-[0 1 1], 'Max error is %0.5e\n', max(normEA));
		
		for j = 1:length(aN)
			dt = t(aN(j)+1)-t(aN(j));
			tAdd(j) = (t(aN(j)+1)+t(aN(j)))/2;
			% Find the x state at the beginning and end of the segment
			x1 = x(28*(aN(j)-1)+1: 28*(aN(j)-1)+7);
			x3 = x(28*(aN(j)-1)+8: 28*(aN(j)-1)+14);
			x5 = x(28*(aN(j)-1)+15: 28*(aN(j)-1)+21);
			x7 = x(28*(aN(j)-1)+22: 28*(aN(j)-1)+28);
			% Find the u at the beginning and end of the segment. Remember that we
			% don't have control defined at the end point, so it is assumed to be
			% the same as before, just for consistency within the code.
			uA = x(nState+4*(aN(j)-1)+1: nState+4*(aN(j)-1)+4);
			
			T = uA(1)*ones(1, 4);
			angleRotated = dt*[tauRatio(1), tauRatio(3), tauRatio(5), tauRatio(7)] + ...
				t(aN(j))*ones(1,4);
			alpha = atan2(uA(3),uA(2))*ones(1, 4) - angleRotated;
			if imag(alpha)
				error('alpha imag')
			end
			beta = asin(uA(4))*ones(1, 4);
			if imag(beta)
				error('beta image')
			end
			Y = [x1, x3, x5, x7];
			Ydot = getDerivCSI(mu, Y, T, alpha, beta, ispND, g0ND);
			
			C = [Y, dt/2*Ydot]*invA;
			xAMat = C*BAdd;
			
			xA = [x1; xAMat(:,1); xAMat(:,2); xAMat(:,3); ...
				xAMat(:,3); xAMat(:, 4); xAMat(:,5); x7];
			
			if length(aN) == 1
				xNew = [x(1:28*(aN(j)-1)); xA; x(28*aN(j)+1: nState)];
				uNew = [x(nState+1:nState+4*(aN(j)-1)+4); uA; ...
					x(nState+4*aN(j)+1:end)];
			else
				if j < length(aN)
					if j == 1
						xNew = [x(1:28*(aN(j)-1)); xA];
						uNew = [x(nState+1: nState+4*(aN(j)-1)+4); uA];
					else
						xNew = [xNew; x(28*aN(j-1)+1:28*(aN(j)-1)); xA];
						uNew = [uNew; x(nState+4*aN(j-1)+1: nState+4*(aN(j)-1)+4); uA];
					end
				else
					xNew = [xNew; x(28*aN(j-1)+1:28*(aN(j)-1)); xA;
						x(28*aN(j)+1:28*s)];
					uNew = [uNew; x(nState+4*aN(j-1)+1: nState+4*(aN(j)-1)+4); uA;
						x(nState+4*aN(j)+1:end)];
				end
			end
		end
		
		t = sort([t; tAdd']);
		StatePhase.timeSegment = t;
		StatePhase.nSegment = length(t)-1;
		nSegment = StatePhase.nSegment;
		StatePhase.state = reshape(xNew, [7, 4*nSegment])';
		StatePhase.control = reshape(uNew, [4, nSegment])';
		doneMeshPhase = false;
		timeVariable = nan(4*nSegment, 1);
		timeDefect = nan(3*nSegment, 1);
		for i = 1:nSegment
			dt = t(i+1) - t(i);
			timeAug = t(i) + dt*tauRatio;
			timeVariable(4*(i-1)+1:4*i) = timeAug(1:2:7);
			timeDefect(3*(i-1)+1:3*i) = timeAug(2:2:6);
		end
		
		StatePhase.timeDefect = timeDefect;
		StatePhase.timeVariable = timeVariable;
		
	else % if aN is empty, or no mesh to add
		doneMeshPhase = true;
	end
end


