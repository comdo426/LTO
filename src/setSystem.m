% setSystem.m* 
%     * Defines system of interest
%     * Defines number of phases
%     * Defines dynamics for each phase

function System = setSystem(nPhase, phaseBody)
% setSystem transforms phase information to system parameters. 
% For nPhase, the input should be a scalar number of phases
% For phaseBody, the input should be a structure with number of bodies

System = cell(nPhase, 1);

for iPhase = 1:nPhase
   
   sizePhasePlanet = size(phaseBody{iPhase});
   mPlanet = sizePhasePlanet(2);
   system{iPhase}.bodies = phaseBody{iPhase};
   
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
   
   %        System{iPhase}.dynamics
   
end % nPhase if loop
    
end

function System2BP = setSystem2BP(P, iPhase)

global MU_SUN
global MU_EARTH
global MU_MOON
global MU_MARS
global SMA_SUNEARTH
global SMA_SUNMARS
global SMA_EARTHMOON
global SMA_MARSDEIMOS
setGlobalVariable

if strcmp(P, 'Earth')
   System2BP.mu = MU_EARTH;
   System2BP.lstar = SMA_EARTHMOON;
   System2BP.tstar = sqrt(SMA_EARTHMOON^3/MU_EARTH);
else if strcmp(P, 'Moon')
      System2BP.mu = MU_MOON;
      System2BP.lstar = SMA_EARTHMOON;
      System2BP.tstar = sqrt(SMA_EARTHMOON^3/MU_MOON);
   else if strcmp(P, 'Sun')
         System2BP.mu = MU_SUN;
         System2BP.lstar = SMA_SUNEARTH;
         System2BP.tstar = sqrt(SMA_SUNEARTH^3/MU_SUN);
      else if strcmp(P, 'Mars')
            System2BP.mu = MU_MARS;
            System2BP.lstar = SMA_MARSDEIMOS;
            System2BP.tstar = sqrt(SMA_MARSDEIMOS^3/MU_MARS);
         else
            fprintf('Phase no. %d: Wrong name of bodies. Only supports S-E-M-Mars', iPhase);
         end
      end
   end
end

end

function SystemCR3BP = setSystemCR3BP(P1, P2, iPhase)

global MU_SUN
global MU_EARTH
global MU_MOON
global MU_MARS
global SMA_SUNEARTH
global SMA_SUNMARS
global SMA_EARTHMOON
global SMA_MARSDEIMOS
setGlobalVariable

if strcmp(P1, 'Earth')
   if strcmp(P2, 'Moon')
      SystemCR3BP.mu = MU_MOON/(MU_EARTH+MU_MOON);
% 		SystemCR3BP.mu = 0.012150586550569; % Override from ATD value
      SystemCR3BP.lstar = SMA_EARTHMOON;
      SystemCR3BP.tstar = sqrt(SMA_EARTHMOON^3/(MU_EARTH+MU_MOON));
		SystemCR3BP.mu1 = MU_EARTH;
		SystemCR3BP.mu2 = MU_MOON;
   else
      fprintf('Phase no. %d: Only supports E-M CR3BP', iPhase);
   end
else if strcmp(P1, 'Sun')
      if strcmp(P2, 'Earth')
         SystemCR3BP.mu = MU_EARTH/(MU_SUN+MU_EARTH);
         SystemCR3BP.lstar = SMA_SUNEARTH;
         SystemCR3BP.tstar = sqrt(SMA_SUNEARTH^3/(MU_SUN+MU_EARTH));
			SystemCR3BP.mu1 = MU_SUN;
			SystemCR3BP.mu2 = MU_EARTH;
      else if strcmp(P2, 'Mars')
            SystemCR3BP.mu = MU_MARS/(MU_SUN+MU_MARS);
            SystemCR3BP.lstar = SMA_SUNMARS;
            SystemCR3BP.tstar = sqrt(SMA_SUNMARS^3/(MU_SUN+MU_MARS));
				SystemCR3BP.mu1 = MU_SUN;
				SystemCR3BP.mu2 = MU_MARS;
         else
            fprintf('Phase no. %d: Only supports S-E, S-M CR3BP', iPhase);
         end
      end
   else
      fprintf('Phase no. %d: P1 should be Earth or Sun', iPhase);
   end
end

      
end

