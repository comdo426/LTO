%TESTSINGLEPHASE - Tests the single-phase collocation with the same
%  dynamics (EM CR3BP)
%
%  Description:
%     Tests the single-phase collocation. This is a down-grade version of the
%     multi-phase collocation, to test the other improvements at the simplest
%     case possible
%
%  Output:
%     Result -	structure of the final state from the solver
%
%  See also: TESTMULTIPHASESAMEDYNAMICS
%	MAT-files required:
%
%  Author: Beom Park
%  Date: 24-Feb-2020; Last revision: 24-Feb-2020

clear; close all;
%% Phase, dynamics information for each phase

nPhase = 1;
phaseBody{1,1} = {'EARTH', 'MOON'};
System = setSystem(nPhase, phaseBody);

%% Initial guess

Setup.transferType = 'periodicOrbits';

% TEST 1: DRO to L4 SPO with atd

% Setup.Phase{1,1}.orbitSource = 'atd';
% Setup.Phase{2,1}.orbitSource = 'atd';
% Setup.Phase{1,1}.orbitName = 'Distant_Retrograde_2D';
% Setup.Phase{1,1}.orbitSelect = {'JC', 2.223};
% Setup.Phase{2,1}.orbitName = 'Short_Period_L4L5';
% Setup.Phase{2,1}.orbitSelect = {'JC', 2.223};

% TEST 2: L2 Halo atd. use Phase 1: s = 17, nRev = 2

Setup.Phase{1,1}.orbitSource = 'atd';
Setup.Phase{1,1}.orbitName = 'Halo_L2';
Setup.Phase{1,1}.orbitSelect = {'JC', 3.015};

% TEST 3: User provided data (DRO.mat, L4SPO.mat). use Phase 1: s = 18, nRev =
% 3, Phase 2: s = 18, nRev = 3 to get the same figures in the example folder
%
% Setup.Phase{1,1}.orbitSource = 'user';
% Setup.Phase{2,1}.orbitSource = 'user';
% Setup.Phase{1,1}.orbitData = load('DRO.mat');
% Setup.Phase{2,1}.orbitData = load('L4SPO.mat');

Setup.Phase{1,1}.s = 17; % Number of segments per rev
Setup.Phase{1,1}.nRev = 3; % Total rev for the orbit stack
% Setup.Phase{2,1}.s = 17;
% Setup.Phase{2,1}.nRev = 2;

LPointVec = [1 2];
Setup.plot = {true, 1, LPointVec};

% plot for the initial guess

InitialGuess = setInitialGuess(System, Setup);
figure(1)
hold on
ini = plot3(InitialGuess{1}.initialConstraint(1), ...
	InitialGuess{1}.initialConstraint(2), ...
	InitialGuess{1}.initialConstraint(3), '^k', 'linewidth', 2.0);
fin = plot3(InitialGuess{end}.finalConstraint(1), ...
	InitialGuess{end}.finalConstraint(2), ...
	InitialGuess{end}.finalConstraint(3), 'vk', 'linewidth', 2.0);
mu = System{1}.parameter.mu;
% 	earthPlot = drawEarth(mu);
moonPlot = drawMoon(mu);
lpPlot = drawLagrangianPoints(mu, Setup.plot{3});

for iPhase = 1:nPhase
	plot3(InitialGuess{iPhase}.state(:,1), ...
		InitialGuess{iPhase}.state(:,2), ...
		InitialGuess{iPhase}.state(:,3), 'b-', 'linewidth', 1.0);
	% TODO - make a function that computes the segment 
	plot3(InitialGuess{iPhase}.state(:,1), ...
		InitialGuess{iPhase}.state(:,2), ...
		InitialGuess{iPhase}.state(:,3), 'k.', 'linewidth', 2.0);
end
legend([ini, fin, moonPlot, lpPlot], ...
	{'Initial', 'Final', 'Moon', 'Li'}, 'fontsize', 13);
xlabel('x(n.d.)', 'fontsize', 13)
ylabel('y(n.d.)', 'fontsize', 13)


%% Spacecraft specs

SC.m0D = 500; % kg
SC.ispD = 2000; % seconds
SC.thrustMaxD = 0.1; % Newton

Spacecraft = setSpacecraft(System, SC);

%% Set optimizer option

Option.integrate = odeset('RelTol', 1e-12, 'AbsTol', 1e-18);
Option.LTO = LToptset( ...
	'FeasibilitySolver', true, ...
	'Optimizer', true, ...
	'MeshRefinement', true);
Option.newton.maxIteration = 200;
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

Option.doneFeasible = false;
Option.doneOptimize = false;
Option.doneMesh = false;
Option.removeMesh = false;
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
% Option.AddCon{1, 1} = {'nlnr' ,'ineq', 'Altitude', 1, 30000};
% Option.AddCon{1, 2} = {'nlnr', 'ineq', 'Altitude', 2, 2000};
%
% Option.AddCon{2, 1} = {'nlnr', 'ineq', 'Altitude', 1, 30000};
% Option.AddCon{2, 2} = {'nlnr', 'ineq', 'Altitude', 2, 7000};

%% Calls LTOMain

Result = LTOMain(System, InitialGuess, Spacecraft, Option);

prompt = 'LTO done, file save?';
fileName = input(prompt);
save(fileName);



