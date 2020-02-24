function [F, dF] = fsolveConstraintEphemerisCoast(x, JDFix, iFix, System, Body, N)
%FSOLVECONSTRAINTEPHEMEIRSCOAST - computes the F and dF matrix for fsolve/newtonRaphson
%
%  Syntax:
%     [F, dF] = FSOLVECONSTRAINTEPHEMERISCOAST(x, JDFix, System, Body)
%
%  Description:
%     Computes F and dF matrix for fsolve/newtonRaphson. Holds for the Ephemeris
%     Coast model
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

nTime = N.time;
nSegment = N.segment;
nState = N.state;

% 6*nSegment: state continuity between segments, nSegment: time continuity
% between segments, 1: JDFix constraints

mState = 6*nSegment;
mTime = nSegment;
F = nan(mState+mTime+1, 1);
dF = zeros(mState+mTime+1, nState+nSegment+nTime);

% time of integration between
tInt = x(nState+1:nState+nSegment);
% time at the beginning/end of segments
t = x(nState+nSegment+1:nState+nSegment+nTime);
% characteristic time
tstar = System.tstar;

for j = 1:nSegment
	
	Y0 = x(6*(j-1)+1:6*j);
	ts = [t(j), tInt(j)+t(j)];
	[~, Y] = ode113(@(t,y) ephemeris(t,y,JDFix, System, Body), ts, Y0, opts);
	h = 1e-8;
	phi = nan(6,6);
	for jj = 1:6
		hvec = zeros(6,1);
		hvec(jj) = h;
		[~, Yp] = ode113(@(t,y) ephemeris(t,y,JDFix, System, Body), ts, Y0+1i*hvec, opts);
% 		[~, Ym] = ode113(@(t,y) ephemeris(t,y,JDFix, System, Body), ts, Y0-1i*hvec, opts);
% 		phi(:,jj) = imag(Yp(end,:)'-Ym(end,:)')/(2*h);
		phi(:,jj) = -imag(Yp(end,:)'-Y(end,:)')/h;
	end
	YEndDeriv = ephemeris(ts(2), Y(end,:)', JDFix, System, Body);
	F(6*(j-1)+1:6*j) = x(6*j+1:6*j+6)-Y(end,:)';
	dF(6*(j-1)+1:6*j, 6*(j-1)+1:6*j) = -phi;
	dF(6*(j-1)+1:6*j, 6*j+1:6*j+6) = eye(6);
	dF(6*(j-1)+1:6*j, nState+j) = -YEndDeriv;
	
	[~, YpTime] = ode113(@(t,y) ephemeris(t,y,JDFix, System, Body), ts+h, Y0, opts);
	[~, YmTime] = ode113(@(t,y) ephemeris(t,y,JDFix, System, Body), ts-h, Y0, opts);
	dF(6*(j-1)+1:6*j, nState+nSegment+j) = -(YpTime(end,:)'-YmTime(end,:)')/(2*h);
	
	F(mState+j) = t(j+1) - (t(j)+tInt(j));
	dF(mState+j, nState+j) = -1;
	dF(mState+j, nState+nSegment+j) = -1;
	dF(mState+j, nState+nSegment+j+1) = 1;
end
% 
F(mState+mTime+1) = t(iFix);
dF(mState+mTime+1, nState+nSegment+iFix) = 1;
	
save('TEST')
end