function [F, dF, dFspr] = fsolveConstraintChariotSEOnly(x, System, State, ...
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
%			state - 7 component for each node, 28 components for each segment
%			control - 4 components(T, ux, uy, uz) for each segment
%			slack variables - 1 component for each segment(Thrust), 2 components
%			for each node(altitude constraint)
%
%  Outputs:
%     F - constraint column vector
%		dF - dF/dx matrix
%
%  See also: FMINCONCONSTRAINT
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 01-Mar-2020

nPhase = length(System);

invA = Collocation.invA;
B = Collocation.B;
D = Collocation.D;


%% number of continuity constraint

mContinuityTTL = 0;
isIniCon = ~isempty(State{1}.initialConstraint);
if isIniCon
	stateIniCon = State{1}.initialConstraint;
	nIniCon = length(stateIniCon);
	mContinuityTTL = mContinuityTTL + nIniCon;
end
% Override isFinCon
isFinCon = 1;
if isFinCon
	nFinCon = 6;
	stateFinCon = State{1}.finalConstraint;
	mContinuityTTL = mContinuityTTL + nFinCon;
end
% isFinCon = ~isempty(State{end}.finalConstraint);
% if isFinCon
% 	stateFinCon = State{end}.finalConstraint;
% 	nFinCon = length(stateFinCon);
% 	mContinuityTTL = mContinuityTTL + nFinCon;
% end
mContinuityTTL = mContinuityTTL + 7*(nPhase-1);

phaseIni = zeros(7, nPhase);
phaseFin = zeros(7, nPhase);

%% number of parameter, number of constraints for each phase

[nx, ~] = getnx(State);
[mc, ~] = getmc(System, State, Spacecraft, Option);

%% dimension of the F, dF matrix

m = sum(mc) + mContinuityTTL; % row number of constraint
n = length(x); % column number of variables
if n ~= sum(nx)
	error('Number of parameters wrong')
end

F = nan(m, 1);
dF = zeros(m, n);

tstar = nan(nPhase,1);

%% iPhase loop

for iPhase = 1:nPhase
	
	%% setting numbers for each phase
	
	isCR3BP = strcmp(System{iPhase}.dynamics(1), 'CR3BP');
	is2BP = strcmp(System{iPhase}.dynamics(1), '2BP');
	mu = System{iPhase}.parameter.mu;
	lstar = System{iPhase}.parameter.lstar;
	tstar(iPhase) = System{iPhase}.parameter.tstar;
	tau0 = State{iPhase}.t0;
	
	[t, nSegment, nNode, nState, nControl, mDefect, ~] = ...
		getPhaseStateInfo(State{iPhase});
	[thrustMaxND, ispND, g0ND] = getSpacecraftInfo(Spacecraft{iPhase});
	% nb, mb: Number of parameters/constraints Before the i-th phase
	if iPhase == 1
		nb = 0;
		mb = 0;
	else
		nb = sum(nx(1:iPhase-1));
		mb = sum(mc(1:iPhase-1));
	end
	
	Y = reshape(x(1+nb:nState+nb), [7,4,nSegment]);
	% TODO - check that this does not break after the first segment;
	uColumn = reshape(x(nState+1+nb:nState+nControl+nb), [4,1,nSegment]);
	UVar = repmat(uColumn, [1,4,1]);
	
	uxVar = UVar(2, :, :);
	uyVar = UVar(3, :, :);
	
	UDef = repmat(uColumn, [1,3,1]);
	
	uxDef = UDef(2, :, :);
	uyDef = UDef(3, :, :);
	
	tVarMat = reshape(State{iPhase}.timeVariable, [1,4,nSegment]);
	tDefMat = reshape(State{iPhase}.timeDefect, [1,3,nSegment]);
	dtMat = nan(1,1,nSegment);
	for i = 1:nSegment
		dt = t(i+1)-t(i);
		dtMat(:,:,i) = dt;
	end
	
	UVarRot = UVar;
	UDefRot = UDef;
	
	UVarRot(2:3, :, :) = [cos(tVarMat+tau0).*uxVar + sin(tVarMat+tau0).*uyVar;
		-sin(tVarMat+tau0).*uxVar + cos(tVarMat+tau0).*uyVar];
	
	UDefRot(2:3, :, :) = [cos(tDefMat+tau0).*uxDef + sin(tDefMat+tau0).*uyDef;
		-sin(tDefMat+tau0).*uxDef + cos(tDefMat+tau0).*uyDef];
	
	Ydot = getDerivCSIVectorized(mu, Y, UVarRot, ispND, g0ND);
	
	C = mtimesx([Y, dtMat/2.*Ydot], invA);
	CB = mtimesx(C,B);
	CD = mtimesx(C,D);
	CBdot = getDerivCSIVectorized(mu, CB, UDefRot, ispND, g0ND);
	FMat = CD - dtMat/2.*CBdot;
	FDefect = reshape(FMat, [21*nSegment, 1]);
	

	h = sqrt(eps);
	
	dF2DMat = nan(21, 28, nSegment);
	
	for jj = 1:28
		hMat = zeros(7,4,nSegment);
		columnNo = ceil(jj/7);
		rowNo = jj-7*(columnNo-1);
		hMat(rowNo, columnNo, :) = ones(1,1,nSegment);
		Yp = Y + h*1j*hMat;
		Ydotp = getDerivCSIVectorized(mu, Yp, UVarRot, ispND, g0ND);
		Cp = mtimesx([Yp, dtMat/2.*Ydotp], invA);
		CBp = mtimesx(Cp,B);
		CDp = mtimesx(Cp,D);
		CBdotp = getDerivCSIVectorized(mu, CBp, UDefRot, ispND, g0ND);
		FMatp = CDp - dtMat/2.*CBdotp;
		dF2DMat(:,jj,:) = reshape(imag(FMatp-FMat)/h, [21, 1, nSegment]);
	end
	
	dFCell = num2cell(dF2DMat, [1 2]);
	dFdx = blkdiag(dFCell{:});
	
	dF2DMatu = nan(21, 4, nSegment);
	
	for jj = 1:4
		
		huMat = zeros(4,1,nSegment);
		huMat(jj,1,:) = ones(1,1,nSegment);
		huColumn = h*1j*huMat;
		UVarp = UVar + repmat(huColumn, [1,4,1]);
		UDefp = UDef + repmat(huColumn, [1,3,1]);
			
		uxVarp = UVarp(2, :, :);
		uyVarp = UVarp(3, :, :);
		
		uxDefp = UDefp(2, :, :);
		uyDefp = UDefp(3, :, :);
			
		UVarRotp = UVarp;
		UDefRotp = UDefp;

		UVarRotp(2:3, :, :) = [cos(tVarMat+tau0).*uxVarp + sin(tVarMat+tau0).*uyVarp;
			-sin(tVarMat+tau0).*uxVarp + cos(tVarMat+tau0).*uyVarp];

		UDefRotp(2:3, :, :) = [cos(tDefMat+tau0).*uxDefp + sin(tDefMat+tau0).*uyDefp;
			-sin(tDefMat+tau0).*uxDefp + cos(tDefMat+tau0).*uyDefp];
	
		Ydotp = getDerivCSIVectorized(mu, Y, UVarRotp, ispND, g0ND);
		
		Cp = mtimesx([Y, dtMat/2.*Ydotp], invA);
		CBp = mtimesx(Cp,B);
		CDp = mtimesx(Cp,D);
		CBdotp = getDerivCSIVectorized(mu, CBp, UDefRotp, ispND, g0ND);
		FMatp = CDp - dtMat/2.*CBdotp;
		dF2DMatu(:,jj,:) = reshape(imag(FMatp-FMat)/h, [21, 1, nSegment]);
	end

	dFCelldu = num2cell(dF2DMatu, [1 2]);
	dFdu = blkdiag(dFCelldu{:});
	
	dFDefect = [dFdx, dFdu, zeros(mDefect, nx(iPhase)-nState-nControl)];
	
	%% Continuity between segments
	
	if nSegment > 1
		FCont = nan(7*(nSegment-1), 1);
		dFCont = zeros(7*(nSegment-1), nx(iPhase));
	else
		FCont = [];
		dFCont = [];
	end
	
	mCont = 7*(nSegment-1);
	
	for i = 1:nSegment-1
		FCont(7*(i-1)+1:7*i) = Y(:,1,i+1) - Y(:,end,i);
		dFCont(7*(i-1)+1:7*i, 28*(i-1)+22:28*(i-1)+28) = -eye(7);
		dFCont(7*(i-1)+1:7*i, 28*i+1:28*i+7) = eye(7);
	end
	%% Thrust slack constraint
	% slack variable at each segment
	mThrust = nSegment;
	nThrust = nSegment;
	FThrust = zeros(mThrust, 1);
	dFThrust = zeros(mThrust, nx(iPhase));
	
	for i = 1:nSegment
		% switch to angle between 0 ~ 2pi just for clarity
		x(nState + nControl + i + nb) = rem(x(nState + nControl + i + nb), 2*pi);
		lambda = x(nState + nControl + i + nb);
		FThrust(i) = uColumn(1,1,i) - thrustMaxND*sin(lambda)^2;
		dFThrust(i, nState + 4*(i-1)+1) = 1;
		dFThrust(i, nState + nControl + i) = ...
			-2*thrustMaxND*sin(lambda)*cos(lambda);
	end
	
	%% Thrust vector unit magnitude constraint
	
	mThrustUnit = nSegment;
	
	FThrustUnit = zeros(mThrustUnit, 1);
	dFThrustUnit = zeros(mThrustUnit, nx(iPhase));
	
	for i = 1:nSegment
		ThrustVec = x(nState + 4*(i-1)+2 + nb: nState + 4*(i-1)+4 + nb);
		FThrustUnit(i) = norm(ThrustVec)^2 - 1;
		dFThrustUnit(i, nState + 4*(i-1)+2: nState + 4*(i-1)+4) = ...
			2*ThrustVec';
	end
	
	% initial/final state at each phase
	phaseIni(:, iPhase) = Y(:,1,1);
	phaseFin(:, iPhase) = Y(:,end,end);
	
	
	%% Additional constraints
	if ~isempty(Option.AddCon{iPhase, 1})
		
		FAlt = nan(2*nNode, 1);
		dFAlt = zeros(2*nNode, nx(iPhase));
		for iNode = 1:nNode
			
			% Note that some of the below is written only for altitude constraint.
			sigma1 = x(nState + nControl + nThrust + iNode + nb);
			nSigma1 = nNode;
			mSigma1 = nNode;
			sigma2 = x(nState + nControl + nThrust + nSigma1 + iNode +nb);
			
			minAlt1 = Option.AddCon{iPhase,1}{5}/lstar;
			minAlt2 = Option.AddCon{iPhase,2}{5}/lstar;
			pos = x(7*(iNode-1)+1+nb:7*(iNode-1)+3+nb);
			alt1 = norm(pos' - [-mu, 0, 0]);
			alt2 = norm(pos' - [1-mu, 0, 0]);
			
			FAlt(iNode, 1) = alt1^2*sin(sigma1) - minAlt1^2;
			FAlt(mSigma1 + iNode, 1) = alt2^2*sin(sigma2) - minAlt2^2;
			
			dFAlt(iNode, 7*(iNode-1)+1:7*(iNode-1)+3) = ...
				2*(pos' - [-mu, 0, 0])*sin(sigma1);
			dFAlt(iNode, nState+nControl+nThrust+iNode) = ...
				alt1^2*cos(sigma1);
			dFAlt(mSigma1 + iNode, 7*(iNode-1)+1:7*(iNode-1)+3) = ...
				2*(pos' - [1-mu, 0, 0])*sin(sigma2);
			dFAlt(mSigma1 + iNode, nState+nControl+nThrust+nSigma1+iNode) = ...
				alt2^2*cos(sigma2);
			
		end
	else
		FAlt = [];
		dFAlt = [];
	end
	F(1+mb:mc(iPhase)+mb) = [FDefect; FCont; FThrust; FThrustUnit; FAlt];
	dF(1+mb:mc(iPhase)+mb, 1+nb:nx(iPhase)+nb) = ...
		[dFDefect; dFCont; dFThrust; dFThrustUnit; dFAlt];
	
	%% Continuity constraint
	
	mt = sum(mc); % number of constraints over the phases
	
	save('C:\Users\comdo\MATLAB Drive\LTOProject\test\test040820\TESTSEOnly')
% 	error(' ')
	
	% F
	switch iPhase
		case 1
			F(mt+1: mt+nIniCon) = phaseIni(:, 1) - stateIniCon;
			F(mt+nIniCon+1:mt+nIniCon+nFinCon) = phaseFin(1:nFinCon,1) - stateFinCon;
		case 2
			% Interphase continuity between EM-SE frame
			EMEndState = phaseFin(:, 1);
			EMEndStateRot = interPhaseRotation(EMEndState, State{1}, State{2}, ...
				Option.FrameSystem.EM, Option.FrameSystem.SE);
			F(mt+nIniCon+7*(iPhase-2)+1: mt+nIniCon+7*(iPhase-1)) = ...
				phaseIni(:, 2) - EMEndStateRot;
		case 3
			% Interphase continuity between SE-SM frame
			SEEndState = phaseFin(:, 2);
			SEEndStateRot = interPhaseRotation(SEEndState, State{2}, State{3}, ...
				Option.FrameSystem.SESC, Option.FrameSystem.SMSC);
			F(mt+nIniCon+7*(iPhase-2)+1: mt+nIniCon+7*(iPhase-1)) = ...
				phaseIni(:, 3) - SEEndStateRot;
			% Continuity for the final state defined by the sprial down
			stateSpiralDown = getStateSpiralDown(System{3}, Spacecraft{3}, ...
				x(end-1), x(end), Option, Option.FrameSystem.SM);
			F(mt+nIniCon+7*(iPhase-1)+1: mt+nIniCon+7*iPhase) = ...
				stateSpiralDown - phaseFin(:,3);			
	end
			
	
	% dF
	
	dF_EMSE = zeros(7,7);
	dF_SESM = zeros(7,7);
	
	switch iPhase
		case 1
			dF(mt+1: mt+nIniCon, 1:nIniCon) = eye(nIniCon);
			dF(mt+nIniCon+1: mt+nIniCon+nFinCon, (nState-7)+1: (nState-7)+ nFinCon) = ...
				eye(nFinCon);


			% Interphase continuity between EM-SE frame
% 			for ii = 1:7
% 				hVec = zeros(7,1);
% 				hVec(ii) = h;
% 				EMEndStatep = phaseFin(:, 1)+hVec;
% 				EMEndStatem = phaseFin(:, 1)-hVec;
% 				EMEndStateRotp = interPhaseRotation(EMEndStatep, State{1}, State{2}, ...
% 					Option.FrameSystem.EM, Option.FrameSystem.SE);
% 				EMEndStateRotm = interPhaseRotation(EMEndStatem, State{1}, State{2}, ...
% 					Option.FrameSystem.EM, Option.FrameSystem.SE);
% 				dF_EMSE(:, ii) = -(EMEndStateRotp-EMEndStateRotm)/(2*h);
% 			end
% 			dF(mt+nIniCon+1: mt+nIniCon+7, ...
% 				(nState-7)+1:nState) = dF_EMSE;
		case 2
			
			dF(mt+nIniCon+1: mt+nIniCon+7, ...
				1+nb: 7+nb) = eye(7);	
			% Interphase continuity between SE-SM frame
			for ii = 1:7
				hVec = zeros(7,1);
				hVec(ii) = h;
				SEEndStatep = phaseFin(:, 2)+hVec;
				SEEndStatem = phaseFin(:, 2)-hVec;
				SEEndStateRotp = interPhaseRotation(SEEndStatep, State{2}, State{3}, ...
					Option.FrameSystem.SESC, Option.FrameSystem.SMSC);
				SEEndStateRotm = interPhaseRotation(SEEndStatem, State{2}, State{3}, ...
					Option.FrameSystem.SESC, Option.FrameSystem.SMSC);
				dF_SESM(:, ii) = -(SEEndStateRotp-SEEndStateRotm)/(2*h);
			end
			dF(mt+nIniCon+7+1: mt+nIniCon+7+7, ...
				(nState-7)+1+nb: nState+nb) = dF_SESM;
			
		case 3
			
			dF(mt+nIniCon+7*(iPhase-2)+1: mt+nIniCon+7*(iPhase-1), ...
				1+nb: 7+nb) = eye(7);
			
			% Continuity for the final state defined by the sprial down back
			% propagation
			dF(mt+nIniCon+7*(iPhase-1)+1: mt+nIniCon+7*iPhase, ...
				(nState-7)+1+nb: nState+nb) = -eye(7);
			
			stateSpiralDownAopp = getStateSpiralDown(System{3}, Spacecraft{3}, ...
				x(end-1)+h, x(end), Option, Option.FrameSystem.SM);
			stateSpiralDownAopm = getStateSpiralDown(System{3}, Spacecraft{3}, ...
				x(end-1)-h, x(end), Option, Option.FrameSystem.SM);
			stateSpiralDownMfp = getStateSpiralDown(System{3}, Spacecraft{3}, ...
				x(end-1), x(end)+h, Option, Option.FrameSystem.SM);
			stateSpiralDownMfm = getStateSpiralDown(System{3}, Spacecraft{3}, ...
				x(end-1), x(end)-h, Option, Option.FrameSystem.SM);
			
			dF(mt+nIniCon+7*(iPhase-1)+1: mt+nIniCon+7*iPhase, ...
				end-1) = (stateSpiralDownAopp-stateSpiralDownAopm)/(2*h);
			dF(mt+nIniCon+7*(iPhase-1)+1: mt+nIniCon+7*iPhase, ...
				end) = (stateSpiralDownMfp-stateSpiralDownMfm)/(2*h);
	end
	
	
end % iPhase for loop
%

% 			save('TESTChariot')
% 			error('TEST')
	

dFspr = sparse(dF);

end