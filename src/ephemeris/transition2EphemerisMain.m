function [statePeriodicOrbitPlot, TransferEphemeris] = ...
	transition2EphemerisMain(System, State, Spacecraft, Option, EphemOption)
%TRANSITION2EPHEMERISMAIN - transition of the converged solution from LTO to
%	ephemeris, main function
%
%  Syntax:
%     output = TRANSITION2EPHEMERISMAIN(System, State, Spacecraft, ...
%		ephemOption)
%
%  Description:
%     transition of the converged solution from LTO to ephemeris. This function
%     summons following three functions
%
%  Inputs:
%		Option - structure that contains following structures
%			ephemOption - structure that contains following information
%				JD0 - initial epoch of the transfer (Julian date, seconds)
%				nRev - number of revolutions stacked for ephemeris conversion of the
%				periodic orbits
%					[1,2] = size(nRev)
%				s - number of segments for each orbit
%					[1,2] = size(s)
%				tOffset - time offset for each orbit to help with ephemeris convegence
%					[1,2] = size(tOffset)
%			newton - newtonRaphson settings
%				maxIteration
%				fcnTolerance
%			fsolve - fsolve settings
%			FrameSystem - CR3BP system, only applicable to CR3BP
%			Body - structure that lists the gravitational parameter and the ID of
%			the bodies participating in the gravity
%
%  Outputs:
%     statePeriodicOR
%
%  See also: TESTTRANSITION2EPHEMERIS
%
%   Author: Beom Park
%   Date: 17-Feb-2020; Last revision: 23-Feb-2020

%% Periodic orbits to ephemeris

cprintf(-[0 0.5 0.5], 'Converging periodic orbits\n');
[statePeriodicOrbitPlot, JDFix, stateFix] = ...
	getEphemerisBoundaryConst(System, State, Option, EphemOption);

%% Transfer arc to ephemeris

cprintf(-[0 0.5 0.5], 'Converging transfer trajectory\n');
TransferEphemeris = getEphemerisTransfer(System, State, Spacecraft, Option,...
	EphemOption, JDFix, stateFix);

end
