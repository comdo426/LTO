function dFstr = ipoptJacobianStructureEphemeris(auxdata)
%IPOPTJACOBIANEPHEMERIS - computes the Jacobian structure of Ephemeris
%Collocation
%
%  Syntax:
%     dFspr = FSOLVECONSTRAINT(x, System, State, ...
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
%   Date: 03-Mar-2020; Last revision: 03-Mar-2020

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
n = sum(nxOpt); % column number of variables
if n ~= sum(nxOpt)
	error('Number of parameters wrong')
end

dF = zeros(m, n);

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
	
	dFDefect = zeros(mDefect, nxOpt(iPhase));
	for i = 1:nSegment
		dFDefect(21*(i-1)+1 :21*(i-1)+21, ...
			28*(i-1)+1:28*(i-1)+28) = ones(21,28);
		dFDefect(21*(i-1)+1 :21*(i-1)+21, ...
			nState+4*(i-1)+1:nState+4*(i-1)+4) = ones(21,4);
	end
	
	%% Continuity between segments
	
	if nSegment > 1
% 		FCont = nan(7*(nSegment-1), 1);
		dFCont = zeros(7*(nSegment-1), nxOpt(iPhase));
	else
		FCont = [];
		dFCont = [];
	end
	
	mCont = 7*(nSegment-1);
	
	for i = 1:nSegment-1
% 		FCont(7*(i-1)+1:7*i) = Y(:,1,i+1) - Y(:,end,i);
		dFCont(7*(i-1)+1:7*i, 28*(i-1)+22:28*(i-1)+28) = eye(7);
		dFCont(7*(i-1)+1:7*i, 28*i+1:28*i+7) = eye(7);
	end
	%% Thrust slack constraint
	% slack variable at each segment
	mThrust = 0;
	nThrust = 0;
% 	FThrust = zeros(mThrust, 1);
	dFThrust = [];
	
		
	%% Thrust vector unit magnitude constraint
	
	mThrustUnit = nSegment;
	
% 	FThrustUnit = zeros(mThrustUnit, 1);
	dFThrustUnit = zeros(mThrustUnit, nxOpt(iPhase));
	
	for i = 1:nSegment
% 		FThrustUnit(i) = norm(ThrustVec)^2 - 1;
		dFThrustUnit(i, nState + 4*(i-1)+2: nState + 4*(i-1)+4) = ...
			[1, 1, 1];
	end
	
	
	%% Additional constraints
	if ~isempty(Option.AddCon{iPhase, 1})
		
		FAlt = nan(2*nNode, 1);
		dFAlt = zeros(2*nNode, nxOpt(iPhase));
		nSigma1 = nNode;
		mSigma1 = nNode;
			
		for iNode = 1:nNode
			
			
			dFAlt(iNode, 7*(iNode-1)+1:7*(iNode-1)+3) = ...
				1;
			dFAlt(iNode, nState+nControl+nThrust+iNode) = ...
				1;
			dFAlt(mSigma1 + iNode, 7*(iNode-1)+1:7*(iNode-1)+3) = ...
				1;
			dFAlt(mSigma1 + iNode, nState+nControl+nThrust+nSigma1+iNode) = ...
				1;
			
		end
	else
		FAlt = [];
		dFAlt = [];
	end
% 	F(1+mb:mcOpt(iPhase)+mb) = [FDefect; FCont; FThrust; FThrustUnit; FAlt];
	dF(1+mb:mcOpt(iPhase)+mb, 1+nb:nxOpt(iPhase)+nb) = ...
		[dFDefect; dFCont; dFThrust; dFThrustUnit; dFAlt];
	
	%% Continuity constraint
	
	mt = sum(mcOpt); % number of constraints over the phases
		
	% dF
	if iPhase == 1
		dF(mt+1: mt+nIniCon, 1:nIniCon) = eye(nIniCon);
		if nPhase > 1
			isSame = isequaln(System{iPhase+1}.dynamics, System{iPhase}.dynamics);
			if isSame
				dF(mt+nIniCon+1: mt+nIniCon+7, (nState-7)+1:nState) = eye(7);
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
				(nState-7)+1+nb: nState+nb) = eye(7);
		end
	end
	
end % iPhase for loop
%

dFstr = sparse(dF);

end