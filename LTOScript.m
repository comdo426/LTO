% LTOScript.m(temporary name)
% Script file to run/test LTOMain
% Dependencies: setSystem, setInitialGuess, setSpacecraft, LTOMain

clear all; close all;



mu1 = 398600.4418;
mu2 = 4902.80107;
mu = mu2/(mu1+mu2);
l_star = 384400;
t_star = sqrt(l_star^3/(mu1+mu2));

L = Lagrange(mu);

%% Get rid of this shit later on

% x0_initial = [1.79949873708526; 0; 0; 0; -1.46727525939563; 0];
% x0_initial = [0.195328385711603; 0.0000000000000; 0; 0; 2.717614643431193 ; 0];
xf_initial = [0.928542597235188; 0.866025403784438; 0; 0.110773126272870; -0.983802039797893; 0];
% xf_initial = [1.79949873708526; 0; 0; 0; -1.46727525939563; 0];

% x0_initial = xf_initial;
x0_initial = [0.195328385711603; 0.0000000000000; 0; 0; 2.717614643431193 ; 0];

P0 = 2*3.135521047051145;
% Pf = P0;
Pf = 6.310163395436344;
% P0 = Pf;
opts = odeset('reltol', 2.3e-14, 'abstol', 1e-19);

% %% One that I used for presentation
n_initial = 3;
n_final = 3;

s_initial = 55;
s_final = 55;

%% One that I used for presentation
% n_initial = 1.5;
% n_final = 1;
% 
% s_initial = 80;
% s_final = 10;

%% New numbers for try(works with 0.2N, sqp)
% n_initial = 5.5;
% n_final = 4;
% s_initial = 30;
% s_final = 20;

tau = [-1; -sqrt(495 +66*sqrt(15))/33; -sqrt(495 -66*sqrt(15))/33; 0; sqrt(495 -66*sqrt(15))/33; sqrt(495 +66*sqrt(15))/33; 1];


t_initial = linspace(0, P0*n_initial, 3*s_initial);
t_s_initial = linspace(0, P0*n_initial, s_initial);
t_final = linspace(0, Pf*n_final, 3*s_final+1);
t_s_final = linspace(0, Pf*n_initial, s_final+1);

x_initial(:,1) = x0_initial;
x_final(:,1) = xf_initial;

for i = 1:length(t_initial)-1
    [T Y] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), [t_initial(i) t_initial(i+1)], x_initial(:,i), opts); 
    x_initial(:,i+1) = Y(end, 1:6)';    
end

% for i = 1:length(t_s_initial)-1
%     dt = t_s_initial(i+1) - t_s_initial(i);
%     [T Y] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), dt/2*[-1 tau(3)], x_initial(:,3*(i-1) +1), opts);
%     x_initial(:,3*(i-1)+2) = Y(end, 1:6)';
%     [T Y] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), dt/2*[tau(3) tau(5)], x_initial(:,3*(i-1) +2), opts);
%     x_initial(:,3*(i-1)+3) = Y(end, 1:6)';
%     [T Y] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), dt/2*[tau(5) tau(7)], x_initial(:,3*(i-1) +3), opts);
%     x_initial(:,3*(i-1)+4) = Y(end, 1:6)';
% end

x_initial = [x_initial; ones(1, 3*s_initial)];

for i = 1:length(t_final)-1
    [T Y] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), [t_final(i) t_final(i+1)], x_final(:,i), opts); 
    x_final(:,i+1) = Y(end, 1:6)';
end

% for i = 1:length(t_s_final)-1
%     dt = t_s_final(i+1) - t_s_final(i);
%     [T Y] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), dt/2*[-1 tau(3)], x_final(:,3*(i-1) +1), opts);
%     x_final(:,3*(i-1)+2) = Y(end, 1:6)';
%     [T Y] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), dt/2*[tau(3) tau(5)], x_final(:,3*(i-1) +2), opts);
%     x_final(:,3*(i-1)+3) = Y(end, 1:6)';
%     [T Y] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), dt/2*[tau(5) tau(7)], x_final(:,3*(i-1) +3), opts);
%     x_final(:,3*(i-1)+4) = Y(end, 1:6)';
% end

x_final = [x_final; ones(1, 3*s_final+1)];
% 
% x_transfer(:, 1) = x_initial(:, end) + (tau(3)+1)/2*(x_final(:,1) - x_initial(:, end));
% x_transfer(:, 2) = x_initial(:, end) + (tau(5)+1)/2*(x_final(:,1) - x_initial(:, end));

s = s_initial + s_final;
n = 3*s;
t_s = [t_s_initial t_s_initial(end)+t_s_initial(end)-t_s_initial(end-1)+t_s_final];

x_DC_initial = zeros(21*s + 7 + 3*s);
Control_initial = [];
for i = 1:s
Control_initial = [Control_initial; 1; 1; 0];
end

x_DC_initial = [reshape([x_initial, x_final], [21*s + 7, 1]); 0.000001*Control_initial];
% x_initial = reshape(x_DC_initial(1:21*s+7), [7, n+1]);w

figure(1)
hold on
plot(x_initial(1, :), x_initial(2, :), 'kx', 'linewidth', 1.5)
plot(x_final(1, :), x_final(2, :), 'kx', 'linewidth', 1.5)
axis equal
figure(2)
hold on
plot(x_initial(4, :), x_initial(5, :), 'bx', 'linewidth', 1.5)
plot(x_final(4, :), x_final(5, :), 'rx', 'linewidth', 1.5)
axis equal




%% Initial guess generation

SystemParameter.mu = mu;
SystemParameter.l_star = l_star;
SystemParameter.t_star = t_star;
SystemParameter.m_0_d = 500; % kg
m_0_d = 500;
SystemParameter.Isp_d = 2000; % Second
SystemParameter.T_max_d = 0.1; % Newton
SystemParameter.P_max = []; % Watt

InitialGuess.x = x_DC_initial;
InitialGuess.x0_initial = x0_initial;
InitialGuess.xf_initial = xf_initial;
InitialGuess.s = s;
InitialGuess.t_s = t_s;

ManifoldsParameter = [];
% ManifoldsParameter
% ManifoldsParameter
% ManifoldsParameter
% ManifoldsParameter

% Thought I had it converged somehow.. but NO!
% load('Optimized_1.mat')
% 
% InitialGuess.x = Result.x_optimized;



%% Optimization Setup

opts = LToptset('EngineType', 'CSI', 'FeasibilitySolver', true, 'Optimizer', false, 'MeshRefinement', true, 'Ephemeris', false, 'NaturalDynamics', 'none');
opts2 = optimoptions('fmincon', 'Algorithm', 'Interior-point', 'Display','iter', 'MaxFunctionEvaluations', 30000000, 'StepTolerance', 1e-20, 'CheckGradient', false, 'SpecifyConstraintGradient', true,  'SpecifyObjectiveGradient', true, 'HessianApproximation', 'bfgs', 'FiniteDifferenceStepSize', 1e-7,'SubproblemAlgorithm', 'cg', 'FiniteDifferenceType', 'central', 'ConstraintTolerance', 1e-13, 'OptimalityTolerance', 1e-7, 'maxiteration', 100000, 'InitBarrierParam', 1e-4);

%% Main function call
Result = LToptimizer_main(SystemParameter, InitialGuess, ManifoldsParameter, opts, opts2);




%% Testing graphs here
% load('Result_meshtest.mat')
