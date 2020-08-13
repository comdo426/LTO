function [F, dF, dFspr] = ...
	fsolveConstraintEphemerisLT(x, JDFix, stateFix, IspND, g0ND, TmaxND, System, Body, N)
%FSOLVECONSTRAINTEPHEMERISCOLLOCATION - computes the F and dF matrix for fsolve/newtonRaphson
%
%  Syntax:
%     [F, dF] = FSOLVECONSTRAINTEPHEMERISCOAST(x, JDFix, System, Body)
%
%  Description:
%     Computes F and dF matrix for fsolve/newtonRaphson. With LT, with
%     collocation
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
%   Date: 18-Feb-2020; Last revision: 18-Feb-2020

opts = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);

nState = N.state;
nTime = N.time;
nSegment = N.segment;
nControl = N.control;
nThrust = N.slack;

% 7*nSegment: state continuity between segments, nSegment: time continuity
% between segments, 2: JDFix constraints at either ends

mState = 7*nSegment;
mTime = nSegment;
mThrust = nSegment;
mThrustUnit = nSegment;
mJD = 2;
mContinuity = 13;
F = nan(mState+mTime+mThrust+mThrustUnit+mJD+mContinuity, 1);
dF = zeros(mState+mTime+mThrust+mThrustUnit+mJD+mContinuity, nState+nSegment+nTime+nControl+nThrust);

% time of integration between
tInt = x(nState+1:nState+nSegment);
% time at the beginning/end of segments
t = x(nState+nSegment+1:nState+nSegment+nTime);
% characteristic time
tstar = System.tstar;



initialJDFix = JDFix.initial;
finalJDFix = JDFix.final;

initialState = stateFix.initial;
finalState = stateFix.final;

for j = 1:nSegment
	
	Y0 = x(7*(j-1)+1:7*j);
	ts = [t(j), tInt(j)+t(j)];
	u = x(nState+nTime+nSegment+4*(j-1)+1:nState+nTime+nSegment+4*j);
	lambda = x(nState+nTime+nSegment+nControl+j);
	
	[~, Y] = ode113(@(t,y) ...
		ephemerisLT(t, y, u, initialJDFix, IspND, g0ND, System, Body), ...
		ts, Y0, opts);
	h = 1e-8;
	
	% numerical partial for the state
	phi = nan(7,7);
	for jj = 1:7
		hvec = zeros(7,1);
		hvec(jj) = h;
		[~, Yp] = ode113(@(t,y) ...
			ephemerisLT(t, y, u, initialJDFix, IspND, g0ND, System, Body), ...
			ts, Y0+hvec, opts);
		[~, Ym] = ode113(@(t,y) ...
			ephemerisLT(t, y, u, initialJDFix, IspND, g0ND, System, Body), ...
			ts, Y0-hvec, opts);
		phi(:,jj) = (Yp(end,:)'-Ym(end,:)')/(2*h);
	end
	
	% partial for the time of integration
	YEndDeriv = ...
		ephemerisLT(ts(2), Y(end,:)', u, initialJDFix, IspND, g0ND, System, Body);
	
	% numerical partial for the control	
	phiControl = nan(7,4);
	for jj = 1:4
		hvec = zeros(4,1);
		hvec(jj) = h;
		[~, Yp] = ode113(@(t,y) ...
			ephemerisLT(t, y, u+hvec, initialJDFix, IspND, g0ND, System, Body), ...
			ts, Y0, opts);
		[~, Ym] = ode113(@(t,y) ...
			ephemerisLT(t, y, u-hvec, initialJDFix, IspND, g0ND, System, Body), ...
			ts, Y0, opts);
		phiControl(:,jj) = (Yp(end,:)'-Ym(end,:)')/(2*h);
	end	
	
	F(7*(j-1)+1:7*j) = x(7*j+1:7*j+7)-Y(end,:)';
	dF(7*(j-1)+1:7*j, 7*(j-1)+1:7*j) = -phi;
	dF(7*(j-1)+1:7*j, 7*j+1:7*j+7) = eye(7);
	dF(7*(j-1)+1:7*j, nState+j) = -YEndDeriv;
	dF(7*(j-1)+1:7*j, nState+nSegment+nTime+4*(j-1)+1:nState+nSegment+nTime+4*j) = -phiControl;
	
	[~, YpTime] = ode113(@(t,y) ...
		ephemerisLT(t, y, u, initialJDFix, IspND, g0ND, System, Body), ...
		ts+h, Y0, opts);
	[~, YmTime] = ode113(@(t,y) ...
		ephemerisLT(t, y, u, initialJDFix, IspND, g0ND, System, Body), ...
		ts-h, Y0, opts);
	dF(7*(j-1)+1:7*j, nState+nSegment+j) = ...
		-(YpTime(end,:)'-YmTime(end,:)')/(2*h);
	
	F(mState+j) = t(j+1) - (t(j)+tInt(j));
	dF(mState+j, nState+j) = -1;
	dF(mState+j, nState+nSegment+j) = -1;
	dF(mState+j, nState+nSegment+j+1) = 1;
	
	F(mState+mTime+j) = u(1) - TmaxND*sin(lambda)^2;
	dF(mState+mTime+j, nState+nSegment+nTime+4*(j-1)+1) = 1;
	dF(mState+mTime+j, nState+nSegment+nTime+nControl+j) = -2*TmaxND*sin(lambda)*cos(lambda);
	
	F(mState+mTime+mThrust+j) = norm(u(2:4))^2-1;
	dF(mState+mTime+mThrust+j, nState+nSegment+nTime+4*(j-1)+2:nState+nSegment+nTime+4*(j-1)+4) = 2*u(2:4)';
	
end
% 
F(mState+mTime+mThrust+mThrustUnit+1) = initialJDFix*60*60*24/tstar + t(1) - ...
	initialJDFix*60*60*24/tstar; % initial JD
F(mState+mTime+mThrust+mThrustUnit+2) = initialJDFix*60*60*24/tstar + t(end) - ...
	finalJDFix*60*60*24/tstar;
dF(mState+mTime+mThrust+mThrustUnit+1, nState+nSegment+1) = 1;
dF(mState+mTime+mThrust+mThrustUnit+2, nState+nSegment+nTime) = 1;

F(mState+mTime+mThrust+mThrustUnit+mJD+1:mState+mTime+mThrust+mThrustUnit+mJD+7) = ...
	x(1:7) - initialState;
F(mState+mTime+mThrust+mThrustUnit+mJD+8:mState+mTime+mThrust+mThrustUnit+mJD+mContinuity) = ... 
	x(nState-6:nState-1) - finalState;
dF(mState+mTime+mThrust+mThrustUnit+mJD+1:mState+mTime+mThrust+mThrustUnit+mJD+7, 1:7) = eye(7);
dF(mState+mTime+mThrust+mThrustUnit+mJD+8:mState+mTime+mThrust+mThrustUnit+mJD+mContinuity, nState-6:nState-1) = eye(6);

% save('TEST3')
dFspr = sparse(dF);
end