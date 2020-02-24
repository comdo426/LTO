function InitialGuess = setInitialGuess(System, Setup)
%SETINITIALGUESS - sets InitialGuess cell for LTO.
%
%  Syntax:
%     InitialGuess = SETSYSTEM(System, Setup)
%
%  Description:
%     Produces InitialGuess(State) cell for LTO using the System cell and the
%     setup cell
%
%  Inputs:
%     System - system information
%			[nPhase,1] = size(System); cell = class(System)
%		Setup - setup for the initialization
%			[nPhase,1]  = size(Setup); cell = class(phaseBody)
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


% if Setup.transferType ... : If you are using manifolds. Room for improvements.

times = 0; % initialization for distinction from the function "times"
nodes = 0;

nPhase = length(System);

InitialGuess = cell(nPhase,1);

for iPhase = 1:nPhase
	s = Setup.Phase{iPhase,1}.s;
	nRev = Setup.Phase{iPhase,1}.nRev;
	isATD = strcmp(Setup.Phase{iPhase,1}.orbitSource, 'atd');
	isUser = strcmp(Setup.Phase{iPhase,1}.orbitSource, 'user');
	
	if isATD
		load(Setup.Phase{iPhase,1}.orbitName)
		isJC = strcmp(Setup.Phase{iPhase,1}.orbitSelect{1}, 'JC');
		nOrbit = length(times);
		if isJC
			% Choose the orbit with JC that is closest to the input JC
			targetJC = Setup.Phase{iPhase,1}.orbitSelect{2};
			targetJCArray = ones(1,nOrbit)*targetJC;
			JCArray = cell2mat(JC);
			diffJCArray = abs(targetJCArray - JCArray);
			[diffJC, orbitIndex] = min(diffJCArray);
			fprintf('Phase %d, Set JC: %d, Orbit JC: %d, difference: %d\n', ...
				iPhase, targetJC, JCArray(orbitIndex), diffJC);
			stateInitial = nodes{orbitIndex}(1,:); % initial condition
			period = times{orbitIndex}(end); % period
		else
			error('Currently only supports JC orbit selection');
		end
	elseif isUser
		stateInitial = Setup.Phase{iPhase,1}.orbitData.stateInitial;
		period = Setup.Phase{iPhase,1}.orbitData.period;
	else
		error('Currently supports atd or user input orbits')
	end
	
	InitialGuess{iPhase,1} = setInitialGuessPhase(s, nRev, stateInitial, period, ...
		System{iPhase});
	
	if iPhase == 1
		InitialGuess{iPhase, 1}.initialConstraint = ...
			[stateInitial, 1]';
	else
		if iPhase == nPhase
			InitialGuess{iPhase, 1}.finalConstraint = ...
				InitialGuess{iPhase,1}.state(end-6:end-1); % open to change, only constrain 6 states
		end
	end
	
	InitialGuess{iPhase, 1}.slack = [];
	InitialGuess{iPhase, 1}.t0 = 0;
	InitialGuess{iPhase, 1}.period = period;
	

end % iPhase for loop

end

