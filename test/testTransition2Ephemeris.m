%TESTTRANSITION2EPHEMERIS - transition of the feasible/optimized solution to
% ephemeris
%
%  Description:
%     Transitions the feasible/optimized solution from LTO. The solution should
%     be in the standard format, i.e., System/State cells from the LTO. This
%     script tests the functions for transitioning the transfers obtained from
%     low-fidelity models.
%			Option.FrameSystem - define the CR3BP system of interest. Currently
%			only supports CR3BP
%			Option.Body - define the bodies that participate in the gravity. 
%			EphemOption - structure with settings about the ephemeris conversion
%				JD - tixed epoch for the initial/final state constrained
%				nRev, s - number of revolutions, number of segments
%				tOffset - trial & error value for the time offset to match the
%				rotated state similar to the one in CR3BP
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

load('L2Halo2L1Halo4Control.mat');
% load('DRO2L4SPO4Control.mat');
% load('L2Halo2L1Halo.mat');
% load('TeenagerTest.mat')
% load('rp=14264.mat')
% Result = InitialGuess;
% load('DRO2L4SPO.mat');
State = Result;
clearvars -except System State Spacecraft

%% Set Option for the transition

Option.integrate = odeset('RelTol', 1e-12, 'AbsTol', 1e-18);
Option.LTO = LToptset( ...
	'FeasibilitySolver', true, ...
	'Optimizer', true, ...
	'MeshRefinement', true);
Option.newton.maxIteration = 200;
Option.newton.fcnTolerance = 1e-11;
Option.newton.stepTolerance = 1e-17;
% Option.newton = []; % un comment this to use fsolve

Option.fsolve = optimoptions( ...
	'fsolve', 'Display','iter', ...
	'MaxFunctionEvaluations', 30000000, ...
	'StepTolerance', 1e-15, ...
	'CheckGradient', true, ... % Usually turn this on to compare gradient
	'FiniteDifferenceType', 'central', ...
	'FiniteDifferenceStepSize', 1e-7, ...
	'SpecifyObjectiveGradient', true, ...
	'OptimalityTolerance', 1e-10, ...
	'maxiteration', 200);

%% Set the optimizer option

% Structure of the partials
Option.ipoptopt.ipopt.jac_c_constant = 'no';
Option.ipoptopt.ipopt.jac_d_constant = 'no';
Option.ipoptopt.ipopt.hessian_constant = 'no';
Option.ipoptopt.ipopt.hessian_approximation = 'limited-memory';
Option.ipoptopt.ipopt.limited_memory_update_type = 'bfgs';

Option.ipoptopt.ipopt.mu_strategy = 'monotone';

Option.ipoptopt.ipopt.max_iter = 3000;

% Set the convergence criteria

% Option.ipoptopt.ipopt.tol = 1e-8;
% Option.ipoptopt.ipopt.acceptable_tol = 5e-5;
Option.ipoptopt.ipopt.acceptable_iter = 1;

Option.ipoptopt.ipopt.constr_viol_tol = 5e-14;
Option.ipoptopt.ipopt.acceptable_constr_viol_tol = 1e-11;

Option.ipoptopt.ipopt.dual_inf_tol = 5e-7;
Option.ipoptopt.ipopt.acceptable_dual_inf_tol = 1e-6;

Option.ipoptopt.ipopt.compl_inf_tol = 1e-3;
Option.ipoptopt.ipopt.acceptable_compl_inf_tol = 1e-2;

Option.ipoptopt.ipopt.mu_init = 1e-1;
% How close can the variables be to the bounds at the beginning
Option.ipoptopt.ipopt.bound_frac = 1e-12;
% How relaxed the bounds are going to be during the optimization
Option.ipoptopt.ipopt.bound_relax_factor = 1e-10;
Option.ipoptopt.ipopt.nlp_scaling_method = 'gradient-based';
% Do you want the outcome to stay within the original bounds?
Option.ipoptopt.ipopt.honor_original_bounds = 'yes';
% Check derivative
% Option.ipoptopt.ipopt.derivative_test = 'first-order';

%% Additional options

% Set the boundaries for the position and velocity
Option.posBoundary = 50000; % km;
Option.velBoundary = 0.5; % km/s;

Option.doneFeasible = false;
Option.doneOptimize = false;
Option.doneMesh = false;
Option.removeMesh = false;

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
% EM.centralBody = 'EARTH';
EM.centralBody = 'MOON';

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
Option.plotNo = [1,1,2,1,5];

Option.isMultipleShooting = false;
%% Additional constraints

Option.AddCon{1, 1} = [];
Option.AddCon{2, 1} = [];


%% EphemOption structure definition

% Change the tOffset to move the boundary constraints to be located at certain
% points in the rotating frame
EphemOption{1,1}.JD = jday(2025, 1, 1, 0, 0, 0)- ...
	(State{2}.timeSegment(end))*tstar/60/60/24;
EphemOption{1,1}.nRev = 30;
EphemOption{1,1}.s = 20;
EphemOption{1,1}.tOffset = 0;

EphemOption{2,1}.JD = jday(2025, 1, 1, 0, 0, 0);
EphemOption{2,1}.nRev = 30;
EphemOption{2,1}.s = 20;
EphemOption{2,1}.tOffset = 0;


%% Calls the main function 

[statePeriodicOrbitPlot, TransferEphemeris] = ...
	transition2EphemerisMain(System, State, Spacecraft, Option, EphemOption);

figure(1)
axis equal
hold on
initial = statePeriodicOrbitPlot.initial;
final = statePeriodicOrbitPlot.final;
initialFix = statePeriodicOrbitPlot.initialFix;
finalFix = statePeriodicOrbitPlot.finalFix;
plot3(initial(:,1), initial(:,2), initial(:,3), 'k');
plot3(final(:,1), final(:,2), final(:,3), 'k');
plot3(initialFix(1,1), initialFix(1,2), initialFix(1,3), 'r^', 'linewidth', 2.0);
plot3(finalFix(1,1), finalFix(1,2), finalFix(1,3), 'rv', 'linewidth', 2.0);

drawMoon(mu);
drawEarth(mu);

c = [0.79, 0.24, 0.87];

transfer = TransferEphemeris.stateConvergedRotPlot;
plot3(transfer(:,1), transfer(:,2), transfer(:,3), 'color', c, 'linewidth', 1.5)

c = [0.47,0.67,0.19];

transfer = TransferEphemeris.stateConvergedRotPlot;
plot3(transfer(:,1), transfer(:,2), transfer(:,3), 'color', c, 'linewidth', 1.5)


		plot3(stateRot(:,1), stateRot(:,2), stateRot(:,3), 'color', c, 'linewidth', 1.5)

