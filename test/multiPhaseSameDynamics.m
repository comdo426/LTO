%MULTIPHASESAMEDYNAMICS - Tests the multi-phase collocation with the same
%  dynamics (EM CR3BP)
%
%  Description:
%     Tests the multi-phase collocation with the same dynamics between the
%     phases. In this particular script, it is set to be EM CR3BP. User can
%     choose either from the existing MAT-file as an initial guess, or choose
%     orbits from the ATD periodic orbits. Note that ATD periodic orbits have
%     limited resolution on Jacobi Constant level, and may result in poor
%     convergence behavior when the catalog fails to deliver the orbits close to
%     the Jacobi Constant.
%
%  Output:
%     Result -	structure of the final state from the solver
%
%  See also: SINGLEPHASE.m,  MULTIPHASEDIFFDYNAMICS.m
%	MAT-files required:
%
%  Author: Beom Park
%  Date: 01-Feb-2020; Last revision: 16-Feb-2020

clear; close all;
%% Phase, dynamics information for each phase

nPhase = 2;
phaseBody{1,1} = {'EARTH', 'MOON'};
phaseBody{2,1} = {'EARTH', 'MOON'};
System = setSystem(nPhase, phaseBody);

%% Initial guess 

Setup.transferType = 'periodicOrbits';

% TEST 1: DRO to L4 SPO

Setup.Phase{1,1}.orbitSource = 'atd';
Setup.Phase{2,1}.orbitSource = 'atd';
% Setup.Phase{1,1}.orbitName = 'Distant_Retrograde_2D';
% Setup.Phase{1,1}.orbitSelect = {'JC', 2.223};
% Setup.Phase{2,1}.orbitName = 'Short_Period_L4L5';
% Setup.Phase{2,1}.orbitSelect = {'nJC', 2.223};

% TEST 2: L2 Halo to L1 Halo

% Setup.Phase{1,1}.orbitName = 'Halo_L2';
% Setup.Phase{1,1}.orbitSelect = {'JC', 3.015};
% Setup.Phase{2,1}.orbitName = 'Halo_L1';
% Setup.Phase{2,1}.orbitSelect = {'JC', 3.015};

% TEST 3: User provided data (L1Halo.mat, L2Halo.mat)

Setup.Phase{1,1}.orbitSource = 'user';
Setup.Phase{2,1}.orbitSource = 'user';
Setup.Phase{1,1}.orbitData = load('L2Halo.mat');
Setup.Phase{2,1}.orbitData = load('L1Halo.mat');

Setup.Phase{1,1}.s = 30; % Number of segments per rev
Setup.Phase{1,1}.nRev = 2; % Total rev for the orbit stack
Setup.Phase{2,1}.s = 30;
Setup.Phase{2,1}.nRev = 2;

LPointVec = [1 2];
Setup.plot = {true, 1, LPointVec};

InitialGuess = setInitialGuess(System, Setup);

% plot for the initial guess
if Setup.plot{1}
	for iPhase = 1:nPhase
		[stateMat, stateSegmentMat, ~] = getStateControlMat(InitialGuess{iPhase});
		figure(Setup.plot{2})
		commonAxisSetting;
		CR3BPAxisSetting;
		plot3(stateSegmentMat(:,1), stateSegmentMat(:,2), stateSegmentMat(:,3), ...
			'k.', 'linewidth', 2.0)
		plot3(stateMat(:,1), stateMat(:,2), stateMat(:,3), 'b-', 'linewidth', 1.0);
	end	
	mu = System{iPhase}.parameter.mu;
% 	earthPlot = drawEarth(mu);
	moonPlot = drawMoon(mu);
	lpPlot = drawLagrangianPoints(mu, Setup.plot{3});
end

%% Spacecraft specs

SC.m0D = 500; % kg
SC.ispD = 2000; % seconds
SC.thrustMaxD = 0.1; % Newton

Spacecraft = setSpacecraft(System, SC);

%% Set optimizer option

Option.integrate = odeset('RelTol', 1e-12, 'AbsTol', 1e-18);
Option.LTO = LToptset( ...
	'FeasibilitySolver', false, ...
	'Optimizer', true, ...
	'MeshRefinement', true);
Option.newton.maxIteration = 1;
Option.newton.fcnTolerance = 1e-12;
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

Option.fmincon = optimoptions( ...
	'fmincon', 'Algorithm', 'Interior-point', ...
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
Option.meshTolerance = 1e-12;

feasPlotNumber = [[11:1:20];[21:1:30]];
optPlotNumber = [[31:1:40];[41:1:50]];
Option.plot.feasible = {true, feasPlotNumber, LPointVec};
Option.plot.optimize = {true, optPlotNumber, LPointVec};

%% Set additional constraints
%
Option.AddCon{1, 1} = [];
Option.AddCon{2, 1} = [];
%
Option.AddCon{1, 1} = {'nlnr' ,'ineq', 'Altitude', 1, 30000};
Option.AddCon{1, 2} = {'nlnr', 'ineq', 'Altitude', 2, 7000};
%
% Option.AddCon{2, 1} = {'nlnr', 'ineq', 'Altitude', 1, 30000};
% Option.AddCon{2, 2} = {'nlnr', 'ineq', 'Altitude', 2, 7000};

%% Calls LTOMain

Result = LTOMain(System, InitialGuess, Spacecraft, Option);


%% Calls Ephemeris?


