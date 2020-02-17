function System = setSystem(nPhase, phaseBody)
%SETSYSTEM - sets System cell for LTO.
%
%  Syntax:
%     System = SETSYSTEM(nPhase, phaseBody)
%
%  Description:
%     Converts phase information into actual system parameters.
%
%  Inputs:
%     nPhase - Number of phase
%			[1,1] = size(nPhase); double = class(target)
%		phaseBody - Bodies participating in each pahse
%			[nPhase,1]  = size(phaseBody); cell = class(phaseBody)
%
%  Outputs:
%     System - cell with system parameters in it.
%			[nPhase,1] = size(System); cell = class(System)
%			System.dynamics - cell with dynamics information
%			System.parameter - struct with system parameters
%				mu - gravitational parameter(2BP, km^2/sec^3)
%					- mass ratio(CR3BP, n.d.)
%				lstar - characteristic length(km)
%				tstar - characteristic time(s)
%				mu1 - gravitational parameter of 1st primary(CR3BP, km^2/sec^3)
%				mu2 - gravitational parameter of 2nd primary(CR3BP, km^2/sec^3)
%
%  Examples:
%     nPhase = 2;
%		phaseBody{1,1} = {'Earth', 'Moon'};
%		phaseBody{2,1} = {'Earth', 'Moon'};
%		System = setSystem(nPhase, phaseBody);
%
%  Other m-files required: SETGLOBALVARIABLE
%  Subfunctions: SETSYSTEM2BP, SETSYSTEMCR3BP
%
%   Author: Beom Park
%   Date: 01-Feb-2020; Last revision: 16-Feb-2020

System = cell(nPhase, 1);

for iPhase = 1:nPhase
	
	sizePhasePlanet = size(phaseBody{iPhase});
	mPlanet = sizePhasePlanet(2);
	
	switch mPlanet
		case 1 % 2BP
			P = phaseBody{iPhase}{1};
			System{iPhase}.dynamics = {'2BP', P};
			System{iPhase}.parameter = setSystem2BP(P, iPhase);
		case 2 % CR3BP
			P1 = phaseBody{iPhase}{1};
			P2 = phaseBody{iPhase}{2};
			System{iPhase}.dynamics = {'CR3BP', P1, P2};
			System{iPhase}.parameter = setSystemCR3BP(P1, P2, iPhase);
		case 3 % BCR3BP
			fprint('Phase no. %d: Only supports up to 3 bodies now', iPhase);
		otherwise
			fprintf('Phase no. %d: Wrong number of bodies. Only supports up to 4 bodies', iPhase);
	end % mPlanet switch loop
	
end % nPhase if loop

end

function System2BP = setSystem2BP(P, iPhase)
%SETSYSTEM2BP - sets System cell for conic(2BP)

global MU_SUN
global MU_EARTH
global MU_MOON
global MU_MARS
global SMA_SUNEARTH
global SMA_SUNMARS
global SMA_EARTHMOON
global SMA_MARSDEIMOS
setGlobalVariable

if strcmp(P, 'EARTH')
	System2BP.mu = MU_EARTH;
	System2BP.lstar = SMA_EARTHMOON;
	System2BP.tstar = sqrt(SMA_EARTHMOON^3/MU_EARTH);
elseif strcmp(P, 'MOON')
	System2BP.mu = MU_MOON;
	System2BP.lstar = SMA_EARTHMOON;
	System2BP.tstar = sqrt(SMA_EARTHMOON^3/MU_MOON);
elseif strcmp(P, 'SUN')
	System2BP.mu = MU_SUN;
	System2BP.lstar = SMA_SUNEARTH;
	System2BP.tstar = sqrt(SMA_SUNEARTH^3/MU_SUN);
elseif strcmp(P, 'MARS')
	System2BP.mu = MU_MARS;
	System2BP.lstar = SMA_MARSDEIMOS;
	System2BP.tstar = sqrt(SMA_MARSDEIMOS^3/MU_MARS);
else
	fprintf('Phase no. %d: Wrong name of bodies. Only supports S-E-M-Mars', iPhase);
end
end

function SystemCR3BP = setSystemCR3BP(P1, P2, iPhase)
%SETSYSTEMCR3BP - sets System cell for CR3BP

global MU_SUN
global MU_EARTH
global MU_MOON
global MU_MARS
global SMA_SUNEARTH
global SMA_SUNMARS
global SMA_EARTHMOON
global SMA_MARSDEIMOS
setGlobalVariable

if strcmp(P1, 'EARTH')
	if strcmp(P2, 'MOON')
		SystemCR3BP.mu = MU_MOON/(MU_EARTH+MU_MOON);
		SystemCR3BP.lstar = SMA_EARTHMOON;
		SystemCR3BP.tstar = sqrt(SMA_EARTHMOON^3/(MU_EARTH+MU_MOON));
		SystemCR3BP.mu1 = MU_EARTH;
		SystemCR3BP.mu2 = MU_MOON;
	else
		fprintf('Phase no. %d: Only supports E-M CR3BP', iPhase);
	end
elseif strcmp(P1, 'SUN')
	if strcmp(P2, 'EARTH')
		SystemCR3BP.mu = MU_EARTH/(MU_SUN+MU_EARTH);
		SystemCR3BP.lstar = SMA_SUNEARTH;
		SystemCR3BP.tstar = sqrt(SMA_SUNEARTH^3/(MU_SUN+MU_EARTH));
		SystemCR3BP.mu1 = MU_SUN;
		SystemCR3BP.mu2 = MU_EARTH;
	elseif strcmp(P2, 'MARS')
		SystemCR3BP.mu = MU_MARS/(MU_SUN+MU_MARS);
		SystemCR3BP.lstar = SMA_SUNMARS;
		SystemCR3BP.tstar = sqrt(SMA_SUNMARS^3/(MU_SUN+MU_MARS));
		SystemCR3BP.mu1 = MU_SUN;
		SystemCR3BP.mu2 = MU_MARS;
	else
		fprintf('Phase no. %d: Only supports S-E, S-M CR3BP', iPhase);
	end
else
	fprintf('Phase no. %d: P1 should be Earth or Sun', iPhase);
end


end

