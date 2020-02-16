%% Phase, dynamics information for each phase

clear all; close all;

nPhase = 4;
phaseBody{1,1} = {'Earth', 'Moon'};
phaseBody{2,1} = {'Earth', 'Moon'};
phaseBody{3,1} = {'Earth', 'Moon'};
phaseBody{4,1} = {'Earth', 'Moon'};
phaseBody{5,1} = {'Earth', 'Moon'};
phaseBody{6,1} = {'Earth', 'Moon'};
phaseBody{7,1} = {'Earth', 'Moon'};
phaseBody{8,1} = {'Earth', 'Moon'};
System = setSystem(nPhase, phaseBody);

%% Initial guess

Setup.transferType = 'periodicOrbits';
% InitialGuessSetup.transferType = 'manifolds';
% 
Setup.Phase{8,1}.orbitName = 'Distant_Retrograde_2D';
Setup.Phase{8,1}.orbitSelect = {'JC', 2.2059};
Setup.Phase{8,1}.s = 30; % Number of segments for each orbit
Setup.Phase{8,1}.nRev = 1; % Total revolution for the thing .. you know

Setup.Phase{7,1}.orbitName = 'Distant_Retrograde_2D';
Setup.Phase{7,1}.orbitSelect = {'JC', 2.4};
Setup.Phase{7,1}.s = 30; % Number of segments for each orbit
Setup.Phase{7,1}.nRev = 1; % Total revolution for the thing .. you know
% 
Setup.Phase{6,1}.orbitName = 'Distant_Retrograde_2D';
Setup.Phase{6,1}.orbitSelect = {'JC', 2.58};
Setup.Phase{6,1}.s = 30; % Number of segments for each orbit
Setup.Phase{6,1}.nRev = 1; % Total revolution for the thing .. you know
% 
Setup.Phase{5,1}.orbitName = 'Distant_Retrograde_2D';
Setup.Phase{5,1}.orbitSelect = {'JC', 2.73};
Setup.Phase{5,1}.s = 20; % Number of segments for each orbit
Setup.Phase{5,1}.nRev = 1; % Total revolution for the thing .. you know
% 
Setup.Phase{4,1}.orbitName = 'Distant_Retrograde_2D';
Setup.Phase{4,1}.orbitSelect = {'JC', 2.79};
Setup.Phase{4,1}.s = 25; % Number of segments for each orbit
Setup.Phase{4,1}.nRev = 1; % Total revolution for the thing .. you know
% 
Setup.Phase{3,1}.orbitName = 'Distant_Retrograde_2D';
Setup.Phase{3,1}.orbitSelect = {'JC', 2.87};
Setup.Phase{3,1}.s = 25; % Number of segments for each orbit
Setup.Phase{3,1}.nRev = 1; % Total revolution for the thing .. you know
%
Setup.Phase{2,1}.orbitName = 'Distant_Retrograde_2D';
Setup.Phase{2,1}.orbitSelect = {'JC', 2.9};
Setup.Phase{2,1}.s = 25; % Number of segments for each orbit
Setup.Phase{2,1}.nRev = 1; % Total revolution for the thing .. you know
%
Setup.Phase{1,1}.orbitName = 'Lyapunov_L1';
Setup.Phase{1,1}.orbitSelect = {'JC', 3.0012};
Setup.Phase{1,1}.s = 25; % Number of segments for each orbit
Setup.Phase{1,1}.nRev = 1; % Total revolution for the thing .. you know

% Setup.Phase{7,1}.orbitName = 'Distant_Retrograde_2D';
% Setup.Phase{7,1}.orbitSelect = {'JC', 2.2059};
% Setup.Phase{7,1}.s = 1; % Number of segments for each orbit
% Setup.Phase{7,1}.nRev = 1; % Total revolution for the thing .. you know
% 
% Setup.Phase{6,1}.orbitName = 'Distant_Retrograde_2D';
% Setup.Phase{6,1}.orbitSelect = {'JC', 2.35};
% Setup.Phase{6,1}.s = 1; % Number of segments for each orbit
% Setup.Phase{6,1}.nRev = 1; % Total revolution for the thing .. you know
% % 
% Setup.Phase{5,1}.orbitName = 'Distant_Retrograde_2D';
% Setup.Phase{5,1}.orbitSelect = {'JC', 2.5};
% Setup.Phase{5,1}.s = 1; % Number of segments for each orbit
% Setup.Phase{5,1}.nRev = 1; % Total revolution for the thing .. you know
% % 
% Setup.Phase{4,1}.orbitName = 'Distant_Retrograde_2D';
% Setup.Phase{4,1}.orbitSelect = {'JC', 2.63};
% Setup.Phase{4,1}.s = 1; % Number of segments for each orbit
% Setup.Phase{4,1}.nRev = 1; % Total revolution for the thing .. you know
% % 
% Setup.Phase{3,1}.orbitName = 'Distant_Retrograde_2D';
% Setup.Phase{3,1}.orbitSelect = {'JC', 2.7};
% Setup.Phase{3,1}.s = 1; % Number of segments for each orbit
% Setup.Phase{3,1}.nRev = 1; % Total revolution for the thing .. you know
% 
% Setup.Phase{2,1}.orbitName = 'Distant_Retrograde_2D';
% Setup.Phase{2,1}.orbitSelect = {'JC', 2.83};
% Setup.Phase{2,1}.s = 1; % Number of segments for each orbit
% Setup.Phase{2,1}.nRev = 1; % Total revolution for the thing .. you know
% %
% Setup.Phase{1,1}.orbitName = 'Distant_Retrograde_2D';
% Setup.Phase{1,1}.orbitSelect = {'JC', 2.9};
% Setup.Phase{1,1}.s = 1; % Number of segments for each orbit
% Setup.Phase{1,1}.nRev = 1; % Total revolution for the thing .. you know
%
LPointVec = [1 2];
Setup.plot = {true, 1, LPointVec};


InitialGuess = setInitialGuess(System, Setup);

%% Spacecraft specs

SC.m0D = 1000; % kg
SC.ispD = 2000; % seconds
SC.thrustMaxD = 0.2; % Newton

Spacecraft = setSpacecraft(System, SC);

%% Set optimizer option

Option.integrate = odeset('RelTol', 1e-12, 'AbsTol', 1e-18);
Option.LTO = LToptset('FeasibilitySolver', true, 'Optimizer', true, ...
	'MeshRefinement', true);
Option.newton.maxIteration = 200;
Option.newton.fcnTolerance = 1e-12;
Option.newton.stepTolerance = 1e-17;
% Option.newton = []; % un comment this to use fsolve

Option.fsolve = optimoptions('fsolve', 'Display','iter', ...
	'MaxFunctionEvaluations', 30000000, ...
	'StepTolerance', 1e-15, ...
	'CheckGradient', true, ...
	'FiniteDifferenceType', 'central', ...
	'FiniteDifferenceStepSize', 1e-7, ...
	'SpecifyObjectiveGradient', true, ...
	'OptimalityTolerance', 1e-10, ...
	'maxiteration', 200);

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
	'OptimalityTolerance', 1e-6, ...
	'maxiteration', 1000000, ...
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
% 
% Option.AddCon{1, 1} = [];
% Option.AddCon{2, 1} = [];
% Option.AddCon{3, 1} = [];
% Option.AddCon{4, 1} = [];
% Option.AddCon{5, 1} = [];
% Option.AddCon{6, 1} = [];
% Option.AddCon{7, 1} = [];
% Option.AddCon{8, 1} = [];

% % 
Option.AddCon{1, 1} = {'nlnr' ,'ineq', 'Altitude', 1, 30000};
Option.AddCon{1, 2} = {'nlnr', 'ineq', 'Altitude', 2, 7000};

Option.AddCon{2, 1} = {'nlnr', 'ineq', 'Altitude', 1, 30000};
Option.AddCon{2, 2} = {'nlnr', 'ineq', 'Altitude', 2, 7000};
Option.AddCon{3, 1} = {'nlnr' ,'ineq', 'Altitude', 1, 30000};
Option.AddCon{3, 2} = {'nlnr', 'ineq', 'Altitude', 2, 7000};
% 
Option.AddCon{4, 1} = {'nlnr', 'ineq', 'Altitude', 1, 30000};
Option.AddCon{4, 2} = {'nlnr', 'ineq', 'Altitude', 2, 7000};
Option.AddCon{5, 1} = {'nlnr' ,'ineq', 'Altitude', 1, 30000};
Option.AddCon{5, 2} = {'nlnr', 'ineq', 'Altitude', 2, 7000};

Option.AddCon{6, 1} = {'nlnr', 'ineq', 'Altitude', 1, 30000};
Option.AddCon{6, 2} = {'nlnr', 'ineq', 'Altitude', 2, 7000};
Option.AddCon{7, 1} = {'nlnr' ,'ineq', 'Altitude', 1, 30000};
Option.AddCon{7, 2} = {'nlnr', 'ineq', 'Altitude', 2, 7000};

Option.AddCon{8, 1} = {'nlnr', 'ineq', 'Altitude', 1, 30000};
Option.AddCon{8, 2} = {'nlnr', 'ineq', 'Altitude', 2, 7000};



%% Calls LTOMain

Result = LTOMain(System, InitialGuess, Spacecraft, Option);


