%TESTTRANSITION2EPHEMERIS - transition of the feasible/optimized solution to
% ephemeris
%
%  Description:
%     Transitions the feasible/optimized solution from LTO. The solution should
%     be in the standard format, i.e., System/State cells from the LTO. This
%     script tests the functions for transitioning the transfers obtained from
%     low-fidelity models.
%
%  Output:
%		TBD
%
%  MAT-files required: 'L2Halo2L1Halo.mat', 'DRO2L4SPO.mat'
%
%  See also: TESTMULTIPHASESAMEDYNAMICS
%
%  Author: Beom Park
%  Date: 17-Feb-2020; Last revision: 23-Feb-2020

clear; close all;

%% Locating the spice toolki, spk, lsk files

if ispc == 1
	addpath('C:\Users\comdo\Desktop\Research\Research_codes\Ephemeris\mice\lib');
	addpath('C:\Users\comdo\Desktop\Research\Research_codes\Ephemeris\mice\src\mice');
	addpath('C:\Users\comdo\Desktop\Research\Research_codes\Ephemeris\mice\data');
	cspice_furnsh({'C:\Users\comdo\Desktop\Research\Research_codes\Ephemeris\EphmData\naif0010.tls.pc','C:\Users\comdo\Desktop\Research\Research_codes\Ephemeris\EphmData\de421.bsp'});
else
	addpath('/Users/beompark/Documents/MATLAB/mice/lib');
	addpath('/Users/beompark/Documents/MATLAB/mice/src/mice');
	addpath('/Users/beompark/Documents/MATLAB/mice/data');
	cspice_furnsh({'/Users/beompark/Documents/MATLAB/Research_codes/Ephemeris/EphmData/naif0010.tls.pc','/Users/beompark/Documents/MATLAB/Research_codes/Ephemeris/EphmData/de421.bsp'});
end

%% Load the converged solution from LTO

load('L2Halo2L1Halo.mat');
% load('DRO2L4SPO.mat');
State = Result;
clearvars -except System State Spacecraft

%% Set Option for the transition

Option.integrate = odeset('RelTol', 1e-12, 'AbsTol', 1e-18);
Option.newton.maxIteration = 200;
Option.newton.fcnTolerance = 1e-10;
Option.newton.stepTolerance = 1e-17;
% Option.newton = []; % uncomment this to use fsolve

Option.fsolve = optimoptions( ...
	'fsolve', 'Display','iter', ...
	'MaxFunctionEvaluations', 30000000, ...
	'StepTolerance', 1e-15, ...
	'CheckGradient', true, ... % Usually turn this on to compare gradient
	'FiniteDifferenceType', 'central', ...
	'FiniteDifferenceStepSize', 1e-8, ...
	'SpecifyObjectiveGradient', true, ...
	'OptimalityTolerance', 1e-10, ...
	'maxiteration', 200);

%% EM structure definition

mu = System{1}.parameter.mu;
lstar = System{1}.parameter.lstar;
tstar = System{1}.parameter.tstar;
mu1 = System{1}.parameter.mu1;
mu2 = System{1}.parameter.mu2;

EM.P1 = 'EARTH';
EM.P2 = 'MOON';
EM.mu = mu;
EM.lstar = lstar;
EM.tstar = tstar;
EM.mu1 = mu1;
EM.mu2 = mu2;
EM.frame = 'ECLIPJ2000';
% Change this to 'MOON' when applicable
EM.centralBody = 'EARTH';
% EM.centralBody = 'MOON';

global MU_SUN
setGlobalVariable;

% Comment/uncomment following lines to make it EM/SEM ephemeris
Body.GM = [mu1, mu2, MU_SUN];
Body.ID = {'EARTH', 'MOON', 'SUN'};
% Body.GM = [mu1, mu2];
% Body.ID = {'EARTH', 'MOON'};
Option.FrameSystem = EM;
Option.Body = Body;
% plotNo, in the order of 
Option.plotNo = [1,1,2,1,5]


%% EphemOption structure definition

% Change the tOffset to move the boundary constraints to be located at certain
% points in the rotating frame
EphemOption{1,1}.JD = jday(2025, 1, 1, 0, 0, 0)- ...
	(State{1}.timeSegment(end)+State{2}.timeSegment(end))*tstar/60/60/24;
EphemOption{1,1}.nRev = 10;
EphemOption{1,1}.s = 30;
EphemOption{1,1}.tOffset = -0.1;

EphemOption{2,1}.JD = jday(2025, 1, 1, 0, 0, 0);
EphemOption{2,1}.nRev = 10;
EphemOption{2,1}.s = 30;
EphemOption{2,1}.tOffset = 0.15;


%% Calls the main function 

[statePeriodicOrbitPlot, TransferEphemeris] = ...
	transition2EphemerisMain(System, State, Spacecraft, Option, EphemOption);




