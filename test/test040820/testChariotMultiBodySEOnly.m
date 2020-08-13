%TESTCHARIOTMULTIBODY - tests the convergence of initial guess for Chariot in
%the medium-fidelity model
%
%  Description: tests the convergence of initial guess for CHARIOT in a
%  medium-fidelity model. This should later be corrected/optimized in ephemeris
%  model
%
%  Output:
%     Result -	structure of the final state from the solver
%
%  See also: TESTMULTIPHASESAMEDYNAMICS
%	MAT-files required:
%
%  Author: Beom Park
%  Date: 04-Apr-2020; Last revision: 04-Apr-2020

clear; close all;
%% Phase, dynamics information for each phase

nPhase = 1;
phaseBody{1,1} = {'SUN', 'EARTH'};
System = setSystem(nPhase, phaseBody);

%% Initial guess 

% load the initial guess file: Note that generating the initial guess is done in
% a separate file for now. This may be integrated later. 

load('InitialGuessChariotSEOnly.mat')


% plot for the initial guess

LPointVec = [];

%% Spacecraft specs

SC.m0D = 180; % kg
SC.ispD = 3000; % seconds
SC.thrustMaxD = 0.06; % Newton

Spacecraft = setSpacecraft(System, SC);

%% Frame System

[EMFrame, SEFrame, SEFrameSC, SMFrame, SMFrameSC] = setChariotFrame;

Option.FrameSystem.EM = EMFrame;
Option.FrameSystem.SE = SEFrame;
Option.FrameSystem.SESC = SEFrameSC;
Option.FrameSystem.SM = SMFrame;
Option.FrameSystem.SMSC = SMFrameSC;

%% Set the feasibility solver option

Option.integrate = odeset('RelTol', 1e-12, 'AbsTol', 1e-18);
JCTarget = 3.000174;
global MU_MARS
a = 23000;
v = sqrt(MU_MARS/a); % dimensional speed from the circular orbit
% Mars-2B arrival condition assuming that argument of periapsis = 0
stateArrivalD = [a; 0; 0; 0; v; 0];
Option.SpiralDown.state = stateArrivalD;
Option.SpiralDown.aop = 3.592756815963031; % rad
Option.SpiralDown.mf = 0.811000000000000; % final mass
Option.SpiralDown.integrate = odeset(Option.integrate, 'Events', ...
	@(t,y) JCTargetFunc(t,y,System{3}.parameter.mu, JCTarget));

Option.LTO = LToptset( ...
	'FeasibilitySolver', true, ...
	'Optimizer', true, ...
	'MeshRefinement', false);
Option.newton.maxIteration = 100;
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

% linear solver scaling
Option.ipoptopt.ipopt.linear_solver = 'mumps';
Option.ipoptopt.ipopt.mumps_pivtol = 1e-10;
% Option.ipoptopt.ipopt.ma57_automatic_scaling = 'yes';

Option.ipoptopt.ipopt.mu_strategy = 'monotone';

Option.ipoptopt.ipopt.max_iter = 3000;

% Set the convergence criteria

% Option.ipoptopt.ipopt.tol = 1e-8;
% Option.ipoptopt.ipopt.acceptable_tol = 5e-5;
Option.ipoptopt.ipopt.acceptable_iter = 1;

Option.ipoptopt.ipopt.constr_viol_tol = 5e-13;
Option.ipoptopt.ipopt.acceptable_constr_viol_tol = 1e-12;

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
Option.posBoundary = 0.1; % n.d.;
Option.velBoundary = 0.1; % n.d.;
Option.aopBoundary = 0.1; % rad;
Option.mfBoundary = 0.0001; % n.d.

Option.doneFeasible = false;
Option.doneOptimize = false;
Option.doneMesh = false;
Option.removeMesh = false;
% mesh refinemeng criteria
Option.meshTolerance = 1e-10;

feasPlotNumber = [[11:1:20];[21:1:30]];
optPlotNumber = [[31:1:40];[41:1:50]];
Option.plot.feasible = {false, feasPlotNumber, LPointVec};
Option.plot.optimize = {false, optPlotNumber, LPointVec};

%% Set additional constraints
%
Option.AddCon{1, 1} = [];
% Option.AddCon{2, 1} = [];
% Option.AddCon{3, 1} = [];
%

%% Calls LTOMain

Result = LTOMainChariotSEOnly(System, InitialGuess, Spacecraft, Option);

%% Plot Definitions

global MU_SUN
global MU_EARTH
global MU_MOON
global MU_MARS
global SMA_EARTHMOON
global SMA_SUNEARTH
global SMA_SUNMARS

setGlobalVariable;

muEM = MU_MOON/(MU_EARTH+MU_MOON);
lstarEM = SMA_EARTHMOON;
tstarEM = sqrt(lstarEM^3/(MU_EARTH+MU_MOON));

muSE = MU_EARTH/(MU_SUN+MU_EARTH);
lstarSE = SMA_SUNEARTH;
tstarSE = sqrt(lstarSE^3/(MU_SUN+MU_EARTH));

muSM = MU_MARS/(MU_SUN+MU_MARS);
lstarSM = SMA_SUNMARS;
tstarSM = sqrt(lstarSM^3/(MU_SUN+MU_MARS));

EMPlotNo = 1;
ECIPlotNo = 2;
SEPlotNo = 4;
SCIPlotNo = 5;
SMPlotNo = 6;

figure(EMPlotNo)
hold on
grid on
axis equal

figure(ECIPlotNo)
hold on
grid on
axis equal

figure(SEPlotNo)
hold on
grid on
axis equal
LSE = Lagrange(muSE);
plot3(LSE(1), 0, 0, 'rx')
plot3(LSE(2), 0, 0, 'rx')
plot3(1-muSE, 0, 0, 'bo')

figure(SCIPlotNo)
hold on
grid on
axis equal
t = [0:0.01:2*pi];
xMD = cos(t).*lstarSM;
yMD = sin(t).*lstarSM;
xED = cos(t).*lstarSE;
yED = sin(t).*lstarSE;
plot(xMD, yMD, '--k')
plot(xED, yED, '--k')

figure(SMPlotNo)
hold on
grid on
axis equal
xM = cos(t).*1;
yM = sin(t).*1;
xE = cos(t).*lstarSE/lstarSM;
yE = sin(t).*lstarSE/lstarSM;
plot(xM, yM, '--k')
plot(xE, yE, '--k')
plot(-muSM, 0, 'o', 'color', [0.9100    0.4100    0.1700], 'linewidth', 5.0);
xlabel('x(n.d. l^*_{SM})', 'fontsize', 13)
ylabel('y(n.d. l^*_{SM})', 'fontsize', 13)
title('SM rotating frame view', 'fontsize', 14)

figure(SCIPlotNo)
hold on
axis equal
grid on

%% Plot part
save('TEST')
% EM Phase
% 
% figure(EMPlotNo)
% for iSegment = 1:Result{1}.nSegment
% 	isThrust = Result{1}.control(iSegment, 1) > Spacecraft{1}.thrustMaxND*0.05;
% 	if isThrust
% 		plot3(Result{1}.state(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
% 			Result{1}.state(4*(iSegment-1)+1:4*(iSegment-1)+4,2), ...
% 			Result{1}.state(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
% 			'r', 'linewidth', 1.0);
% 	else
% 		plot3(Result{1}.state(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
% 			Result{1}.state(4*(iSegment-1)+1:4*(iSegment-1)+4,2), ...
% 			Result{1}.state(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
% 			'b', 'linewidth', 1.0);
% 	end
% end

% SE Phase

figure(SEPlotNo)
for iSegment = 1:Result{1}.nSegment
	isThrust = Result{1}.control(iSegment, 1) > Spacecraft{1}.thrustMaxND*0.05;
	if isThrust
		plot3(Result{1}.state(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
			Result{1}.state(4*(iSegment-1)+1:4*(iSegment-1)+4,2), ...
			Result{1}.state(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
			'r', 'linewidth', 1.0);
	else
		plot3(Result{1}.state(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
			Result{1}.state(4*(iSegment-1)+1:4*(iSegment-1)+4,2), ...
			Result{1}.state(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
			'b', 'linewidth', 1.0);
	end
end

SESCI = rot2inert(Result{1}.state, Result{1}.timeVariable*tstarSE, 0, SEFrameSC, Result{1}.t0);

figure(SCIPlotNo)
for iSegment = 1:Result{1}.nSegment
	isThrust = Result{1}.control(iSegment, 1) > Spacecraft{1}.thrustMaxND*0.05;
	if isThrust
		plot3(SESCI(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
			SESCI(4*(iSegment-1)+1:4*(iSegment-1)+4,2), ...
			SESCI(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
			'r', 'linewidth', 1.0);
	else
		plot3(SESCI(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
			SESCI(4*(iSegment-1)+1:4*(iSegment-1)+4,2), ...
			SESCI(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
			'b', 'linewidth', 1.0);
	end
end


% SM Phase
% 
% figure(SMPlotNo)
% for iSegment = 1:Result{3}.nSegment
% 	isThrust = Result{3}.control(iSegment, 1) > Spacecraft{3}.thrustMaxND*0.05;
% 	if isThrust
% 		plot3(Result{3}.state(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
% 			Result{3}.state(4*(iSegment-1)+1:4*(iSegment-1)+4,2), ...
% 			Result{3}.state(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
% 			'r', 'linewidth', 1.0);
% 	else
% 		plot3(Result{3}.state(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
% 			Result{3}.state(4*(iSegment-1)+1:4*(iSegment-1)+4,2), ...
% 			Result{3}.state(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
% 			'b', 'linewidth', 1.0);
% 	end
% end
% 
% SMSCI = rot2inert(Result{3}.state, Result{3}.timeVariable*tstarSM, 0, SMFrameSC, Result{3}.t0);
% 
% figure(SCIPlotNo)
% for iSegment = 1:Result{3}.nSegment
% 	isThrust = Result{3}.control(iSegment, 1) > Spacecraft{3}.thrustMaxND*0.05;
% 	if isThrust
% 		plot3(SMSCI(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
% 			SMSCI(4*(iSegment-1)+1:4*(iSegment-1)+4,2), ...
% 			SMSCI(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
% 			'r', 'linewidth', 1.0);
% 	else
% 		plot3(SMSCI(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
% 			SMSCI(4*(iSegment-1)+1:4*(iSegment-1)+4,2), ...
% 			SMSCI(4*(iSegment-1)+1:4*(iSegment-1)+4,1), ...
% 			'b', 'linewidth', 1.0);
% 	end
% end
