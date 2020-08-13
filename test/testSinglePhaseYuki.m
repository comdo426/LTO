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

nPhase = 5;
phaseBody{1,1} = {'EARTH', 'MOON'};
phaseBody{2,1} = {'EARTH', 'MOON'};
phaseBody{3,1} = {'EARTH', 'MOON'};
phaseBody{4,1} = {'EARTH', 'MOON'};
phaseBody{5,1} = {'EARTH', 'MOON'};
System = setSystem(nPhase, phaseBody);
System{1}.parameter.mu = 0.0121506642671889;
System{2}.parameter.mu = 0.0121506642671889;
System{3}.parameter.mu = 0.0121506642671889;
System{4}.parameter.mu = 0.0121506642671889;
System{5}.parameter.mu = 0.0121506642671889;

%% Initial guess

Setup.transferType = 'periodicOrbits';

% TEST 1: DRO to L4 SPO with atd

Setup.Phase{1,1}.orbitSource = 'user';
Setup.Phase{2,1}.orbitSource = 'user';
Setup.Phase{3,1}.orbitSource = 'user';
Setup.Phase{4,1}.orbitSource = 'user';
Setup.Phase{5,1}.orbitSource = 'user';

Setup.Phase{1,1}.orbitData.stateInitial = ...
[1.17068416907643;0;-0.0917764936422443;0;-0.191133784351078;0]';
Setup.Phase{1,1}.orbitData.period = 3.34252831519353;

Setup.Phase{2,1}.orbitData.stateInitial = ...
[1.15801334271845;0;-0.129538449619934;0;-0.210972062069635;0]';
Setup.Phase{2,1}.orbitData.period = 3.25466937068710;

Setup.Phase{3,1}.orbitData.stateInitial = ...
[1.14789580032591;0;-0.150376533060736;0;-0.219571379483693;0]';
Setup.Phase{3,1}.orbitData.period = 3.17659370825530;

Setup.Phase{4,1}.orbitData.stateInitial = ...
[1.13966236457917;0;-0.163884170729207;0;-0.223550264711580;0]';
Setup.Phase{4,1}.orbitData.period = 3.10509451032695;

Setup.Phase{5,1}.orbitData.stateInitial = ...
[1.13266616020310;0;-0.173446411568361;0;-0.225213015739275;0]';
Setup.Phase{5,1}.orbitData.period = 3.03759103899258;

Setup.Phase{1,1}.s = 50; % Number of segments per rev
Setup.Phase{1,1}.nRev = 1; % Total rev for the orbit stack
Setup.Phase{2,1}.s = 50; % Number of segments per rev
Setup.Phase{2,1}.nRev = 1; % Total rev for the orbit stack
Setup.Phase{3,1}.s = 50; % Number of segments per rev
Setup.Phase{3,1}.nRev = 1; % Total rev for the orbit stack
Setup.Phase{4,1}.s = 50; % Number of segments per rev
Setup.Phase{4,1}.nRev = 1; % Total rev for the orbit stack
Setup.Phase{5,1}.s = 50; % Number of segments per rev
Setup.Phase{5,1}.nRev = 1; % Total rev for the orbit stack
% Setup.Phase{2,1}.s = 17;
% Setup.Phase{2,1}.nRev = 2;

LPointVec = [1 2];
Setup.plot = {true, 1, LPointVec};

% plot for the initial guess

InitialGuess = setInitialGuess(System, Setup);
% figure(1)
% hold on
% ini = plot3(InitialGuess{1}.initialConstraint(1), ...
% 	InitialGuess{1}.initialConstraint(2), ...
% 	InitialGuess{1}.initialConstraint(3), '^k', 'linewidth', 2.0);
% fin = plot3(InitialGuess{end}.finalConstraint(1), ...
% 	InitialGuess{end}.finalConstraint(2), ...
% 	InitialGuess{end}.finalConstraint(3), 'vk', 'linewidth', 2.0);
% mu = System{1}.parameter.mu;
% % 	earthPlot = drawEarth(mu);
% moonPlot = drawMoon(mu);
% lpPlot = drawLagrangianPoints(mu, Setup.plot{3});
% 
% for iPhase = 1:nPhase
% 	plot3(InitialGuess{iPhase}.state(:,1), ...
% 		InitialGuess{iPhase}.state(:,2), ...
% 		InitialGuess{iPhase}.state(:,3), 'b-', 'linewidth', 1.0);
% 	% TODO - make a function that computes the segment 
% 	index = getIndexSegment(InitialGuess{iPhase});
% 	plot3(InitialGuess{iPhase}.state(index,1), ...
% 		InitialGuess{iPhase}.state(index,2), ...
% 		InitialGuess{iPhase}.state(index,3), 'k.', 'linewidth', 2.0);
% end
% legend([ini, fin, moonPlot, lpPlot], ...
% 	{'Initial', 'Final', 'Moon', 'Li'}, 'fontsize', 13);
% xlabel('x(n.d.)', 'fontsize', 13)
% ylabel('y(n.d.)', 'fontsize', 13)


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

% Option.fmincon = optimoptions( ...
% 	'fmincon', 'Algorithm', 'Interior-point', ...
% 	'Display','iter', ...
% 	'MaxFunctionEvaluations', 30000000, ...
% 	'StepTolerance', 1e-20, ...
% 	'CheckGradient', false, ...
% 	'SpecifyConstraintGradient', true, ...
% 	'SpecifyObjectiveGradient', true,  ...
% 	'HessianApproximation', 'bfgs', ...
% 	'FiniteDifferenceStepSize', 1e-7, ...
% 	'SubproblemAlgorithm', 'cg', ...
% 	'FiniteDifferenceType', 'central', ...
% 	'ConstraintTolerance', 1e-13, ...
% 	'OptimalityTolerance', 1e-6, ...
% 	'maxiteration', 1000000, ...
% 	'InitBarrierParam', 1e-4);
Option.ipoptopt.ipopt.jac_c_constant = 'no';
Option.ipoptopt.ipopt.jac_d_constant = 'no';

Option.ipoptopt.ipopt.hessian_constant = 'no';
Option.ipoptopt.ipopt.hessian_approximation = 'limited-memory';
Option.ipoptopt.ipopt.limited_memory_update_type = 'bfgs';

Option.ipoptopt.ipopt.mu_strategy = 'monotone';

Option.ipoptopt.ipopt.max_iter = 3000;
Option.ipoptopt.ipopt.tol = 1e-8;
Option.ipoptopt.ipopt.acceptable_tol = 5e-5;

Option.ipoptopt.ipopt.mu_init = 1e-1;
Option.ipoptopt.ipopt.bound_frac = 1e-14;
Option.ipoptopt.ipopt.bound_relax_factor = 1e-10;
Option.ipoptopt.ipopt.nlp_scaling_method = 'gradient-based';
Option.ipoptopt.ipopt.acceptable_constr_viol_tol = 1e-13;

Option.ipoptopt.ipopt.honor_original_bounds = 'yes';
% 
% Option.ipoptopt.ipopt.acceptable_dual_inf_tol = 1e-14;
% Option.ipoptopt.ipopt.derivative_test = 'first-order';

Option.posBoundary = 10000; % km;
Option.velBoundary = 0.5; % km/s;


% Option.ipoptopt.ipopt.derivative_test = 'first-order';





Option.doneFeasible = false;
Option.doneOptimize = false;
Option.doneMesh = false;
Option.removeMesh = false;
Option.meshTolerance = 1e-10;

feasPlotNumber = [[11:1:20];[21:1:30]];
optPlotNumber = [[31:1:40];[41:1:50]];
Option.plot.feasible = {true, feasPlotNumber, LPointVec};
Option.plot.optimize = {true, optPlotNumber, LPointVec};

%% Set additional constraints
%
Option.AddCon{1, 1} = [];
Option.AddCon{2, 1} = [];
Option.AddCon{3, 1} = [];
Option.AddCon{4, 1} = [];
Option.AddCon{5, 1} = [];
%
% Option.AddCon{1, 1} = {'nlnr' ,'ineq', 'Altitude', 1, 30000};
% Option.AddCon{1, 2} = {'nlnr', 'ineq', 'Altitude', 2, 2000};
%
% Option.AddCon{2, 1} = {'nlnr', 'ineq', 'Altitude', 1, 30000};
% Option.AddCon{2, 2} = {'nlnr', 'ineq', 'Altitude', 2, 7000};

%% Calls LTOMain

tic
Result = LTOMain(System, InitialGuess, Spacecraft, Option);
toc

% 
% prompt = 'LTO done, file save?';
% fileName = input(prompt);
% save(fileName);



