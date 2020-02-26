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
%   Date: 01-Feb-2020; Last revision: 24-Feb-2020

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
	uColumn = reshape(x(nState+1+nb:nState+nControl+nb), [3,1,nSegment]);
	UVar = repmat(uColumn, [1,4,1]);
	UDef = repmat(uColumn, [1,3,1]);
	
	tVarMat = reshape(State{iPhase}.timeVariable, [1,4,nSegment]);
	tDefMat = reshape(State{iPhase}.timeDefect, [1,3,nSegment]);
	dtMat = nan(1,1,nSegment);
	for i = 1:nSegment
		dt = t(i+1)-t(i);
		dtMat(:,:,i) = dt;
	end
	
	UVar(2,:,:) = UVar(2,:,:) - tVarMat;
	UDef(2,:,:) = UDef(2,:,:) - tDefMat;
	 
	Ydot = getDerivCSIVectorized(mu, Y, UVar, ispND, g0ND);
	
	C = mtimesx([Y, dtMat/2.*Ydot], invA);
	CB = mtimesx(C,B);
	CD = mtimesx(C,D);
	CBdot = getDerivCSIVectorized(mu, CB, UDef, ispND, g0ND);
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
		Ydotp = getDerivCSIVectorized(mu, Yp, UVar, ispND, g0ND);
		Cp = mtimesx([Yp, dtMat/2.*Ydotp], invA);
		CBp = mtimesx(Cp,B);
		CDp = mtimesx(Cp,D);
		CBdotp = getDerivCSIVectorized(mu, CBp, UDef, ispND, g0ND);
		FMatp = CDp - dtMat/2.*CBdotp;
		dF2DMat(:,jj,:) = reshape(imag(FMatp-FMat)/h, [21, 1, nSegment]);
	end
	
	dFCell = num2cell(dF2DMat, [1 2]);
	dFdx = blkdiag(dFCell{:});
	
	dF2DMatu = nan(21, 3, nSegment);
	
	for jj = 1:3
		huMat = zeros(3,1,nSegment);
		huMat(jj,1,:) = ones(1,1,nSegment);
		uColumnp = uColumn+h*1j*huMat;
		UVarp = repmat(uColumnp, [1,4,1]);
		UDefp = repmat(uColumnp, [1,3,1]);
		
		UVarp(2,:,:) = UVarp(2,:,:) - tVarMat;
		UDefp(2,:,:) = UDefp(2,:,:) - tDefMat;
		
		Ydotp = getDerivCSIVectorized(mu, Y, UVarp, ispND, g0ND);
		
		Cp = mtimesx([Y, dtMat/2.*Ydotp], invA);
		CBp = mtimesx(Cp,B);
		CDp = mtimesx(Cp,D);
		CBdotp = getDerivCSIVectorized(mu, CBp, UDefp, ispND, g0ND);
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
	% 	lambdaVec = x(nState + nControl + 1 + nb : nState + nControl + nThrust + nb);
	%
	% 	ThrustVec = reshape(uColumn(1,1,:), [nSegment, 1]);
	% 	FThrust = ThrustVec - thrustMaxND.*sin(lambdaVec).^2;
	FThrust = zeros(mThrust, 1);
	dFThrust = zeros(mThrust, nx(iPhase));
	
	for i = 1:nSegment
		lambda = x(nState + nControl + i + nb);
		FThrust(i) = uColumn(1,1,i) - thrustMaxND*sin(lambda)^2;
		dFThrust(i, nState + 3*(i-1)+1) = 1;
		dFThrust(i, nState + nControl + i) = ...
			-2*thrustMaxND*sin(lambda)*cos(lambda);
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
	save('TEST')
	F(1+mb:mc(iPhase)+mb) = [FDefect; FCont; FThrust; FAlt];
	dF(1+mb:mc(iPhase)+mb, 1+nb:nx(iPhase)+nb) = [dFDefect; dFCont; dFThrust; dFAlt];
	
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
			% Julian Day to be added from the previous time
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

end