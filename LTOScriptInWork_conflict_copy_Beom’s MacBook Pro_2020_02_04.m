%% Phase, dynamics information for each phase

clear all; close all;

nPhase = 2;
phaseBody{1,1} = {'Earth', 'Moon'};
phaseBody{2,1} = {'Earth', 'Moon'};
System = setSystem(nPhase, phaseBody);

%% Initial guess

Setup.transferType = 'periodicOrbits';
% InitialGuessSetup.transferType = 'manifolds';

% Setup.Phase{1,1}.orbitName = 'Distant_Retrograde_2D';
% Setup.Phase{1,1}.orbitSelect = {'JC', 2.223};
Setup.Phase{1,1}.orbitName = 'Halo_L2';
Setup.Phase{1,1}.orbitSelect = {'JC', 3.015};
Setup.Phase{1,1}.s = 17; % Number of segments for each orbit
Setup.Phase{1,1}.nRev = 3; % Total revolution for the thing .. you know

Setup.Phase{2,1}.orbitName = 'Halo_L1';
Setup.Phase{2,1}.orbitSelect = {'JC', 3.015};
% Setup.Phase{2,1}.orbitName = 'Short_Period_L4L5';
% Setup.Phase{2,1}.orbitSelect = {'JC', 2.223};
Setup.Phase{2,1}.s = 17; % Number of segments for each orbit
Setup.Phase{2,1}.nRev = 3; % Total revolution for the thing .. you know

LPointVec = [1 2 3 4 5];
Setup.plot = {true, 1, LPointVec};

InitialGuess = setInitialGuess(System, Setup);

%% Spacecraft specs

SC.m0D = 500; % kg
SC.ispD = 2000; % seconds
SC.thrustMaxD = 0.1; % Newton

Spacecraft = setSpacecraft(System, SC);

%% Set optimizer option

Option.integrate = odeset('RelTol', 1e-12, 'AbsTol', 1e-18);
Option.LTO = LToptset('FeasibilitySolver', true, 'Optimizer', false, ...
	'MeshRefinement', true);
Option.newton.maxIteration = 200;
Option.newton.fcnTolerance = 1e-12;
Option.newton.stepTolerance = 1e-17;
% Option.newton = []; % un comment this to use fsolve

Option.fsolve = optimoptions('fsolve', 'Display','iter', ...
	'MaxFunctionEvaluations', 30000000, ...
	'StepTolerance', 1e-15, ...
	'CheckGradient', false, ...
	'FiniteDifferenceType', 'central', ...
	'FiniteDifferenceStepSize', 1e-7, ...
	'SpecifyObjectiveGradient', true, ...
	'OptimalityTolerance', 1e-10, ...
	'maxiteration', 300);

Option.fmincon = optimoptions('fmincon', 'Algorithm', 'Interior-point', ...
	'Display','iter', ...
	'MaxFunctionEvaluations', 30000000, ...
	'StepTolerance', 1e-20, ...
	'CheckGradient', false, ...
	'SpecifyConstraintGradient', true, ...
	'SpecifyObjectiveGradient', true,  ...
	'HessianApproximation', 'bfgs', ...
	'FiniteDifferenceStepSize', 1e-7, ...
	'SubproblemAlgorithm', 'cg', ...
	'FiniteDifferenceType', 'central', ...
	'ConstraintTolerance', 1e-13, ...
	'OptimalityTolerance', 1e-7, ...
	'maxiteration', 100000, ...
	'InitBarrierParam', 1e-4);

Option.doneFeasible = 0;
Option.doneOptimize = 0;
Option.doneMesh = 0;
Option.removeMesh = 0;
Option.meshTolerance = 1e-11;

feasPlotNumber = [[11:1:20];[21:1:30]];
optPlotNumber = [[31:1:40];[41:1:50]];
Option.plot.feasible = {true, feasPlotNumber, LPointVec};
Option.plot.optimize = {true, optPlotNumber, LPointVec};

%% Set additional constraints

% Option.AddCon{1, 1} = [];
% Option.AddCon{2, 1} = [];

Option.AddCon{1, 1} = {'nlnr' ,'ineq', 'Altitude', 1, 60000};
Option.AddCon{1, 2} = {'nlnr', 'ineq', 'Altitude', 2, 7000};

Option.AddCon{2, 1} = {'nlnr', 'ineq', 'Altitude', 1, 60000};
Option.AddCon{2, 2} = {'nlnr', 'ineq', 'Altitude', 2, 7000};





%% Calls LTOMain

Result = LTOMain(System, InitialGuess, Spacecraft, Option);

% load('Newton')
% testDefect(mu, x1, x3, x5, x7, u, ispND, g0ND, dt, phi, phiPrime)

% Problem = setProblem(InitialGuess, Spacecraft, Option);


