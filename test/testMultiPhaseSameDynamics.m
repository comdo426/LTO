%TESTMULTIPHASESAMEDYNAMICS - Tests the multi-phase collocation with the same
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
%  See also: TESTTRANSITION2EPHEMERIS
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

% TEST 1: DRO to L4 SPO with atd 25/3
% 
% Setup.Phase{1,1}.orbitSource = 'atd';
% Setup.Phase{2,1}.orbitSource = 'atd';
% Setup.Phase{1,1}.orbitName = 'Distant_Retrograde_2D';
% Setup.Phase{1,1}.orbitSelect = {'JC', 2.223};
% Setup.Phase{2,1}.orbitName = 'Short_Period_L4L5';
% Setup.Phase{2,1}.orbitSelect = {'JC', 2.223};

% TEST 2: L2 Halo to L1 Halo atd. use Phase 1: s = 17, nRev =
% 2, Phase 2: s = 17, nRev = 2 to get the same figures in the example folder

Setup.Phase{1,1}.orbitSource = 'atd';
Setup.Phase{2,1}.orbitSource = 'atd';
Setup.Phase{1,1}.orbitName = 'Halo_L2';
Setup.Phase{1,1}.orbitSelect = {'JC', 3.015180};
Setup.Phase{2,1}.orbitName = 'Halo_L1';
Setup.Phase{2,1}.orbitSelect = {'JC', 3.016348};

% TEST 3: User provided data (DRO.mat, L4SPO.mat). use Phase 1: s = 18, nRev =
% 3, Phase 2: s = 18, nRev = 3 to get the same figures in the example folder

% Setup.Phase{1,1}.orbitSource = 'user';
% Setup.Phase{2,1}.orbitSource = 'user';
% Setup.Phase{1,1}.orbitData = load('DRO.mat');
% Setup.Phase{2,1}.orbitData = load('L4SPO.mat');

% TEST 4: L2 Halo to L1 Halo atd. use Phase 1: s = 17, nRev =
% 2, Phase 2: s = 17, nRev = 2 to get the same figures in the example folder
% %
% Setup.Phase{1,1}.orbitSource = 'atd';
% Setup.Phase{2,1}.orbitSource = 'atd';
% Setup.Phase{1,1}.orbitName = 'Halo_L2';
% Setup.Phase{1,1}.orbitSelect = {'JC', 3.0347};
% Setup.Phase{2,1}.orbitName = 'Distant_Retrograde_2D';
% Setup.Phase{2,1}.orbitSelect = {'JC', 2.9328};


Setup.Phase{1,1}.s = 18; % Number of segments per rev
Setup.Phase{1,1}.nRev = 2.1; % Total rev for the orbit stack
Setup.Phase{2,1}.s = 18;
Setup.Phase{2,1}.nRev = 2.1;

LPointVec = [1 2];
Setup.plot = {true, 1, LPointVec};

% plot for the initial guess

InitialGuess = setInitialGuess(System, Setup);
figure(1)
axis equal
hold on
ini = plot3(InitialGuess{1}.initialConstraint(1), ...
	InitialGuess{1}.initialConstraint(2), ...
	InitialGuess{1}.initialConstraint(3), '^k', 'linewidth', 2.0);
fin = plot3(InitialGuess{2}.finalConstraint(1), ...
	InitialGuess{2}.finalConstraint(2), ...
	InitialGuess{2}.finalConstraint(3), 'vk', 'linewidth', 2.0);

mu = System{1}.parameter.mu;
moonPlot = drawMoon(mu);
lpPlot = drawLagrangianPoints(mu, Setup.plot{3});

for iPhase = 1:nPhase
	plot3(InitialGuess{iPhase}.state(:,1), ...
		InitialGuess{iPhase}.state(:,2), ...
		InitialGuess{iPhase}.state(:,3), 'b-', 'linewidth', 1.0);
	% TODO - make a function that computes the segment 
	index = getIndexSegment(InitialGuess{iPhase});
	plot3(InitialGuess{iPhase}.state(index,1), ...
		InitialGuess{iPhase}.state(index,2), ...
		InitialGuess{iPhase}.state(index,3), 'k.', 'linewidth', 2.0);
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

%% Set the feasibility solver option

Option.integrate = odeset('RelTol', 1e-12, 'AbsTol', 1e-18);
Option.LTO = LToptset( ...
	'FeasibilitySolver', true, ...
	'Optimizer', true, ...
	'MeshRefinement', true);
Option.newton.maxIteration = 100;
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
Option.ipoptopt.ipopt.acceptable_constr_viol_tol = 1e-13;

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
Option.velBoundary = 0.1; % km/s;

Option.doneFeasible = false;
Option.doneOptimize = false;
Option.doneMesh = false;
Option.removeMesh = false;
% mesh refinemeng criteria
Option.meshTolerance = 1e-11;

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
Option.AddCon{1, 2} = {'nlnr', 'ineq', 'Altitude', 2, 10000};
%
Option.AddCon{2, 1} = {'nlnr', 'ineq', 'Altitude', 1, 30000};
Option.AddCon{2, 2} = {'nlnr', 'ineq', 'Altitude', 2, 7000};

%% Calls LTOMain

Result = LTOMain(System, InitialGuess, Spacecraft, Option);

% prompt = 'LTO done, file save?';
% fileName = input(prompt);
% save(fileName);



