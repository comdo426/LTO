function [StatePhase, doneMeshPhase] = ...
	CEPMeshRefineSinglePhase(SystemPhase, StatePhase, SpacecraftPhase, Option)

isCR3BP = strcmp(SystemPhase.dynamics{1}, 'CR3BP');
is2BP = strcmp(SystemPhase.dynamics{1}, '2BP');

t0 = StatePhase.t0;

opts = Option.integrate;
tol = Option.meshTolerance;

mu = SystemPhase.parameter.mu;
[t, s, ~, nState, ~, ~, ~] = getPhaseStateInfo(StatePhase);
[~, ispND, g0ND] = getSpacecraftInfo(SpacecraftPhase);

x = [StatePhase.state; StatePhase.control];

[~, x_s, u] = getStateControlMat(StatePhase);

tau = [-1; -sqrt(495 +66*sqrt(15))/33; -sqrt(495 -66*sqrt(15))/33; 0; ...
	sqrt(495 -66*sqrt(15))/33; sqrt(495 +66*sqrt(15))/33; 1];
tauRatio = (tau - ones(1, 7)*tau(1))/2;
tauAdd = [-1; -1+0.5*(1-sqrt(495 -66*sqrt(15))/33); ...
	-1+0.5*(1+sqrt(495 -66*sqrt(15))/33); 0; ...
	+ 0.5*(1-sqrt(495 -66*sqrt(15))/33);  0.5*(1+sqrt(495 -66*sqrt(15))/33); 1];

[Phi, PhiPrime, PhiMeshAdd] = LGL_7th_coefficient;
if Option.removeMesh == 1
	
	%% Remove
	
	k = 1;
	rN = []; % short word for remove number
	skipSwitch = 0;
	normER = nan(s-1, 1);
	for i = 1:s-1
		if skipSwitch ~= 1
			% First propagation
			[~, state] = ode113(@(t,y) CR3BPLT_CSI_Inertial(t,y,u(i,:),mu, ...
				ispND,g0ND,t0), [t(i), t(i+1)], x_s(i,:), opts);
			% Second propagation
			[~, state] = ode113(@(t,y) CR3BPLT_CSI_Inertial(t,y,u(i,:),mu, ...
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
			xj1 = x(21*(rN(j)-1)+1:21*(rN(j)-1)+7);
			xj2 = x(21*(rN(j)+1)+1:21*(rN(j)+1)+7);
			xR(:, j) = [xj1+(xj2-xj1)*tauRatio(3); xj1+(xj2-xj1)*tauRatio(5)];
			if nRemove == 1
				xNew = [x(1:21*(rN(j)-1)+7); xR; x(21*(rN(j)+1)+nState)];
				uNew = [x(nState+1: nState+3*(rN(j)-1)+3); x(nState+3*nR(j)+4:end)];
			else
				if j < nRemove
					if j == 1
						xNew = [x(1:21*(rN(j)-1)+7); xR];
						uNew = [x(nState+1: nState+3*(rN(j)-1)+3)];
					else
						xNew = [xNew; x(21*(nR(j-1)+1)+1: 21*(nR(j)-1)+7); xR];
						uNew = [uNew; x(nState+3*(nR(j-1)+1)+1: nState+3*(nR(j)-1)+3)];
					end
				else
					xNew = [xNew; x(21*(nR(j-1)+1)+1: 21*(nR(j)-1)+7); xR; ...
						x(21*(nR(j)+1)+1: nState)];
					uNew = [uNew; x(nState+3*(nR(j-1)+1)+1: nState+3*(nR(j)-1)+3); ...
						x(nState+3*nR(j)+4:end)];
				end
			end
		end
		StatePhase.state = xNew;
		StatePhase.control = uNew;
		t(rN(:)+1) = [];
		StatePhase.timeSegment = t;
		StatePhase.nSegment = length(t)-1;
	end
	
	doneMeshPhase = 0;
else
	Option.removeMesh = 0;
	
	%% Add
	k = 1;
	aN = [];
	normEA = nan(s, 1);
	for i = 1:s
		[~, state] = ode113(@(t,y) CR3BPLT_CSI_Inertial(t,y,u(i,:),mu,ispND, ...
			g0ND, t0), [t(i), t(i+1)], x_s(i,:), opts);
		xProp = state(end,:);
		E = xProp - x_s(i+1, :);
		normEA(i) = norm(E);
		if normEA(i) > tol
			aN(k) = i;
			fprintf('Add: Segment no. %d, the error is %0.5e\n', ...
				aN(k), normEA(i));
			k = k+1;
		end
	end
	
	if ~isempty(aN)
		fprintf('Max error is %0.5e\n', max(normEA));
		
		for j = 1:length(aN)
			dt = t(aN(j)+1)-t(aN(j));
			tAdd(j) = (t(aN(j)+1)+t(aN(j)))/2;
			% Find the x state at the beginning and end of the segment
			x1 = x(21*(aN(j)-1)+1: 21*(aN(j)-1)+7);
			x3 = x(21*(aN(j)-1)+8: 21*(aN(j)-1)+14);
			x5 = x(21*(aN(j)-1)+15: 21*(aN(j)-1)+21);
			x7 = x(21*(aN(j)-1)+22: 21*(aN(j)-1)+28);
			% Find the u at the beginning and end of the segment. Remember that we
			% don't have control defined at the end point, so it is assumed to be
			% the same as before, just for consistency within the code.
			uA = x(nState+3*(aN(j)-1)+1: nState+3*(aN(j)-1)+3);
			
			T = uA(1)*ones(1, 4);
			angleRotated = dt*[tauRatio(1), tauRatio(3), tauRatio(5), tauRatio(7)] + ...
				t(aN(j))*ones(1,4);
			alpha = uA(2)*ones(1, 4) - angleRotated;
			beta = uA(3)*ones(1, 4);
			Y = [x1, x3, x5, x7];
			Ydot = getDerivCSI(mu, Y, T, alpha, beta, ispND, g0ND);
			
			B = [Y, dt/2*Ydot];
			
			xAMat = B*PhiMeshAdd';
			xA = [xAMat(:,1); xAMat(:,2); xAMat(:,3); xAMat(:, 4); xAMat(:,5)];
			
			if length(aN) == 1
				xNew = [x(1:21*(aN(j)-1)+7); xA; x(21*aN(j)+1: nState)];
				uNew = [x(nState+1:nState+3*(aN(j)-1)+3); uA; ...
					x(nState+3*aN(j)+1:end)];
			else
				if j < length(aN)
					if j == 1
						xNew = [x(1:21*(aN(j)-1)+7); xA];
						uNew = [x(nState+1: nState+3*(aN(j)-1)+3); uA];
					else
						xNew = [xNew; x(21*aN(j-1)+1:21*(aN(j)-1)+7); xA];
						uNew = [uNew; x(nState+3*aN(j-1)+1: nState+3*(aN(j)-1)+3); uA];
					end
				else
					xNew = [xNew; x(21*aN(j-1)+1:21*(aN(j)-1)+7); xA;
						x(21*aN(j)+1:21*s+7)];
					uNew = [uNew; x(nState+3*aN(j-1)+1: nState+3*(aN(j)-1)+3); uA;
						x(nState+3*aN(j)+1:end)];
				end
			end
		end
		
		StatePhase.state = xNew;
		StatePhase.control = uNew;
		size(t)
		size(tAdd)
		t = sort([t; tAdd']);
		StatePhase.timeSegment = t;
		StatePhase.nSegment = length(t)-1;
		doneMeshPhase = 0;
	else % if aN is empty, or no mesh to add
		doneMeshPhase = 1;
	end
end


