function F = ipoptConstraintEphemeris(x, auxdata)
%IPOPTCONSTRAINT - computes the constraint vector for ipopt
%
%  Syntax:
%     F = ipoptGradient(x, auxdata)
%
%  Description:
%     computes the constraint vector for ipopt.
%
%  Inputs:
%		x - parameter composed of state, control and slack variables
%		auxdata - the structure that contains
%
%  Outputs:
%     F - constraint vector 
%
%	See also: FSOLVECONSTRAINT, IPOPTJACOBIAN
%
%   Author: Beom Park
%   Date: 26-Feb-2020; Last revision: 01-Mar-2020

System = auxdata.System;
State = auxdata.State;
Spacecraft = auxdata.Spacecraft;
Option = auxdata.Option;
Collocation = auxdata.Collocation;

nPhase = length(System);

invA = Collocation.invA;
B = Collocation.B;
D = Collocation.D;

FrameSystem = Option.FrameSystem;
Body = Option.Body;
CB = FrameSystem.centralBody;
iCB = find(contains(Body.ID, CB), 1);
if isempty(iCB)
	error('cannot find central body in the Body structure')
end
% Note that this version only supports single-phase for the frame selection
muND = Body.GM*FrameSystem.tstar^2/FrameSystem.lstar^3;


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

[~, nxOpt] = getnx(State);
[~, mcOpt] = getmc(System, State, Spacecraft, Option);

%% dimension of the F, dF matrix

m = sum(mcOpt) + mContinuityTTL; % row number of constraint
n = length(x); % column number of variables
if n ~= sum(nxOpt)
	error('Number of parameters wrong')
end

F = nan(m, 1);

tstar = nan(nPhase,1);

%% iPhase loop

for iPhase = 1:nPhase
	
	%% setting numbers for each phase
	bodyStateMatVar = State{iPhase}.bodyStateMatVar;
	bodyStateMatDef = State{iPhase}.bodyStateMatDef;
	
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
		nb = sum(nxOpt(1:iPhase-1));
		mb = sum(mcOpt(1:iPhase-1));
	end
	
	Y = reshape(x(1+nb:nState+nb), [7,4,nSegment]);
	% TODO - check that this does not break after the first segment;
	uColumn = reshape(x(nState+1+nb:nState+nControl+nb), [4,1,nSegment]);
	UVar = repmat(uColumn, [1,4,1]);
	
	UDef = repmat(uColumn, [1,3,1]);
	
	tVarMat = reshape(State{iPhase}.timeVariable, [1,4,nSegment]);
	tDefMat = reshape(State{iPhase}.timeDefect, [1,3,nSegment]);
	dtMat = nan(1,1,nSegment);
	for i = 1:nSegment
		dt = t(i+1)-t(i);
		dtMat(:,:,i) = dt;
	end
	
	Ydot = getDerivCSIVectorizedEphemeris(muND, Y, UVar, ispND, g0ND, bodyStateMatVar, iCB);
	
	C = mtimesx([Y, dtMat/2.*Ydot], invA);
	CB = mtimesx(C,B);
	CD = mtimesx(C,D);
	CBdot = getDerivCSIVectorizedEphemeris(muND, CB, UDef, ispND, g0ND, bodyStateMatDef, iCB);
	FMat = CD - dtMat/2.*CBdot;
	FDefect = reshape(FMat, [21*nSegment, 1]);
	
	
	%% Continuity between segments
	if nSegment > 1
		FCont = nan(7*(nSegment-1), 1);
	else
		FCont = [];
	end
	
	for i = 1:nSegment-1
		FCont(7*(i-1)+1:7*i) = Y(:,1,i+1) - Y(:,end,i);
	end
	
	
	%% Thrust slack constraint
	% slack variable at each segment
	
	nThrust = 0;
	FThrust = [];
% 	
% 	for i = 1:nSegment
% 		lambda = x(nState + nControl + i + nb);
% 		FThrust(i) = uColumn(1,1,i) - thrustMaxND*sin(lambda)^2;
% 	end
	
	
	%% Thrust vector unit magnitude constraint
	
	mThrustUnit = nSegment;
	
	FThrustUnit = zeros(mThrustUnit, 1);
	
	for i = 1:nSegment
		ThrustVec = x(nState + 4*(i-1)+2 + nb: nState + 4*(i-1)+4 + nb);
		FThrustUnit(i) = norm(ThrustVec)^2 - 1;
	end
	
	
	% initial/final state at each phase
	phaseIni(:, iPhase) = Y(:,1,1);
	phaseFin(:, iPhase) = Y(:,end,end);
	
	
	%% Additional constraints
	if ~isempty(Option.AddCon{iPhase, 1})
						
		FAlt = nan(2*nNode, 1);
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
						
		end
	else
		FAlt = [];
	end
	
	%% Merge	
	
	F(1+mb:mcOpt(iPhase)+mb) = [FDefect; FCont; FThrust; FThrustUnit; FAlt];
	%% Continuity constraint
	
	mt = sum(mcOpt);
	
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
	
end % iPhase for loop
% 
% save('TEST3')
end