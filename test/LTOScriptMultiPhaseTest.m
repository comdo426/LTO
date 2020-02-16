%% Phase, dynamics information for each phase

clear; close all;

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

nPhase = 2;
phaseBody{1,1} = {'Earth', 'Moon'};
phaseBody{2,1} = {'Sun'};
System = setSystem(nPhase, phaseBody);

%% Initial guess

LPointVec = [1 2];

% InitialGuess = setInitialGuess(System, Setup);
% load('LICInitialGuess.mat');
load('LICInitialGuess25.mat');
% load('LICInitialGuessGradientCheck.mat');
% load('LICInitialGuessWOThrust.mat');
%% Spacecraft specs
SC.m0D = 1000; % kg
SC.ispD = 2000; % seconds
SC.thrustMaxD = 0.1; % Newton

Spacecraft = setSpacecraft(System, SC);

%% Set optimizer option

Option.integrate = odeset('RelTol', 1e-12, 'AbsTol', 1e-18);
Option.LTO = LToptset('FeasibilitySolver', true, 'Optimizer', true, ...
	'MeshRefinement', false);
Option.newton.maxIteration = 30;
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
	'OptimalityTolerance', 1e-5, ...
	'maxiteration', 1000000 , ...
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
Option.AddCon{1, 1} = [];
Option.AddCon{2, 1} = [];
% % 
Option.AddCon{1, 1} = {'nlnr' ,'ineq', 'Altitude', 1, 10000};
Option.AddCon{1, 2} = {'nlnr', 'ineq', 'Altitude', 2, 1000};

% Option.AddCon{2, 1} = {'nlnr', 'ineq', 'Altitude', 1, 30000};
% Option.AddCon{2, 2} = {'nlnr', 'ineq', 'Altitude', 2, 7000};





%% Calls LTOMain

Result = LTOMain(System, InitialGuess, Spacecraft, Option);


%% Calls Ephemeris? 


