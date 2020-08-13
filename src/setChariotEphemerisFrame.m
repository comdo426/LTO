function [EMFrame, EMFrameMC, SEFrame, SEFrameSC, SMFrame, SMFrameSC] = setChariotEphemerisFrame


global MU_SUN
global MU_EARTH
global MU_MOON
global MU_MARS
global SMA_EARTHMOON
global SMA_SUNEARTH
global SMA_SUNMARS

setGlobalVariable;

muEM = MU_MOON/(MU_EARTH+MU_MOON);
lstarEM = SMA_EARTHMOON;
tstarEM = sqrt(lstarEM^3/(MU_EARTH+MU_MOON));

muSE = MU_EARTH/(MU_SUN+MU_EARTH);
lstarSE = SMA_SUNEARTH;
tstarSE = sqrt(lstarSE^3/(MU_SUN+MU_EARTH));

muSM = MU_MARS/(MU_SUN+MU_MARS);
lstarSM = SMA_SUNMARS;
tstarSM = sqrt(lstarSM^3/(MU_SUN+MU_MARS));

EMFrame.P1 = 'EARTH';
EMFrame.P2 = 'MOON';
EMFrame.mu1 = MU_EARTH;
EMFrame.mu2 = MU_MOON;
EMFrame.centralBody = 'EARTH';
EMFrame.frame = 'ECLIPJ2000';
EMFrame.lstar = lstarEM;
EMFrame.tstar = tstarEM;
EMFrame.mu = muEM;

EMFrameMC = EMFrame;
EMFrameMC.centralBody = 'MOON';

SEFrame.P1 = 'SUN';
SEFrame.P2 = 'EARTH';
SEFrame.mu1 = MU_SUN;
SEFrame.mu2 = MU_EARTH;
SEFrame.centralBody = 'EARTH';
SEFrame.frame = 'ECLIPJ2000';
SEFrame.lstar = lstarSE;
SEFrame.tstar = tstarSE;
SEFrame.mu = muSE;

SEFrameSC = SEFrame;
SEFrameSC.centralBody = 'SUN';

SMFrame.P1 = 'SUN';
SMFrame.P2 = 'MARS';
SMFrame.mu1 = MU_SUN;
SMFrame.mu2 = MU_MARS;
SMFrame.centralBody = 'MARS';
SMFrame.frame = 'ECLIPJ2000';
SMFrame.lstar = lstarSM;
SMFrame.tstar = tstarSM;
SMFrame.mu = muSM;

SMFrameSC = SMFrame;
SMFrameSC.centralBody = 'SUN';


end