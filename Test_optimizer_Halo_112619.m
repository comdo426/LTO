clear all; close all;

mu1 = 398600.4418;
mu2 = 4902.80107;
mu = mu2/(mu1+mu2);
l_star = 384400;
t_star = sqrt(l_star^3/(mu1+mu2));

L = Lagrange(mu);

%% Get rid of this shit later on

% Original example

% x0_initial_guess = [318014.12/l_star; 0; 34596/l_star; 0; 0.2146/l_star*t_star; 0];
% % xf_initial_guess = [449363.6/l_star; 0; 36997/l_star; 0; -0.1994/l_star*t_star; 0];
% xf_initial_guess = [449363.6/l_star; 0; 36998/l_star; 0; -0.1994/l_star*t_star; 0];
% P0_initial_guess = 11.9711*24*3600/t_star;
% Pf_initial_guess = 14.4809*24*3600/t_star;
% 
% % [x0_initial P0_half k F] = Halo_corrector(mu, x0_initial_guess, P0_initial_guess/2);
% % P0 = P0_half*2;
% % [xf_initial Pf_half k F] = Halo_corrector(mu, xf_initial_guess, Pf_initial_guess/2);
% % Pf = Pf_half*2;
% 
% [xf_initial P0_half k F] = Halo_corrector(mu, x0_initial_guess, P0_initial_guess/2);
% Pf = P0_half*2;
% [x0_initial Pf_half k F] = Halo_corrector(mu, xf_initial_guess, Pf_initial_guess/2);
% P0 = Pf_half*2;

% One from Robert's thesis

xf_initial = [0.840645749967222;                   0;   0.159210178951199;                   0;   0.261992933019641;                   0];
Pf = 2*1.349018321829528;

x0_initial =   [1.115001963068551;                   0;   0.190879918445258;                   0;  -0.223386275592327;                   0];
P0 = 2*1.418376500453653;



% P0 = Pf;
opts = odeset('reltol', 2.3e-14, 'abstol', 1e-19);

% %% One that I used for presentation
n_initial = 2;
n_final = 2;

s_initial = 50;
s_final = 50;

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

% load('112719_mesh_examination.mat')
% load('112619_optimized_success_finite.mat');
load('Inertial_test.mat')
InitialGuess.x = Result.x_optimized;
InitialGuess.s = length(Result.t_s_feasible)-1;
InitialGuess.t_s = Result.t_s_feasible;

% for i = 1:InitialGuess.s
%     InitialGuess.x(21*InitialGuess.s + 7 + 3*(i-1) + 2) = InitialGuess.x(21*InitialGuess.s + 7 + 3*(i-1) + 2) - InitialGuess.t_s(i);
% end
%% Optimization Setup

opts = LToptset('EngineType', 'CSI', 'FeasibilitySolver', true, 'Optimizer', true, 'MeshRefinement', true, 'Ephemeris', false, 'NaturalDynamics', 'none');
% opts2 = optimoptions('fmincon', 'Algorithm', 'Interior-point', 'Display','iter', 'MaxFunctionEvaluations', 30000000, 'StepTolerance', 1e-20, 'CheckGradient', false, 'SpecifyConstraintGradient', true,  'SpecifyObjectiveGradient', true, 'HessianApproximation', 'finite-difference', 'FiniteDifferenceStepSize', 1e-7,'SubproblemAlgorithm', 'cg', 'FiniteDifferenceType', 'central', 'ConstraintTolerance', 1e-13, 'OptimalityTolerance', 1e-7, 'maxiteration', 100000, 'InitBarrierParam', 1e-10);
opts2 = optimoptions('fmincon', 'Algorithm', 'Interior-point', 'Display','iter', 'MaxFunctionEvaluations', 30000000, 'StepTolerance', 1e-20, 'CheckGradient', false, 'SpecifyConstraintGradient', true,  'SpecifyObjectiveGradient', true, 'HessianApproximation', 'bfgs', 'FiniteDifferenceStepSize', 1e-7,'SubproblemAlgorithm', 'cg', 'FiniteDifferenceType', 'central', 'ConstraintTolerance', 1e-13, 'OptimalityTolerance', 1e-7, 'maxiteration', 100000, 'InitBarrierParam', 1e-5);
% opts2 = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display','iter', 'MaxFunctionEvaluations', 30000000, 'StepTolerance', 1e-20, 'CheckGradient', false, 'SpecifyConstraintGradient', true,  'SpecifyObjectiveGradient', true, 'FiniteDifferenceType', 'central', 'ConstraintTolerance', 1e-13, 'OptimalityTolerance', 2e-7, 'maxiteration', 100000);


%% Main function call
Result = LToptimizer_main(SystemParameter, InitialGuess, ManifoldsParameter, opts, opts2);



%% plot part
% k = 1;
% 
% x_feasible = reshape(Result.x_optimized(1:21*s +7), [7, n+1]);
% for i = 1:s+1
%     x_feasible_segment(:,i) = Result.x_optimized(21*(i-1) +1: 21*(i-1)+7);
% end
% u_feasible = reshape(Result.x_optimized(21*s + 7 +1:21*s+7+3*s), [3, s]);
% figure(50+k)
% hold on
% for i = 1:s
%     if u_feasible(1,i) < 0.0001
%         color = [0 0 1];
%     else
%         color = [1 0 0];
%     end
%     p_feasible = plot3(x_feasible(1, 3*(i-1)+1:3*i+1), x_feasible(2, 3*(i-1)+1:3*i+1), x_feasible(3, 3*(i-1)+1:3*i+1), 'Color', color, 'linewidth', 1.5);
% end
% p4 = plot3(x_feasible_segment(1, :), x_feasible_segment(2,:), x_feasible_segment(3,:), 'k.', 'linewidth', 1.0);
% axis equal
% plot3(L(1), 0, 0, 'rx', 'linewidth', 1.5);
% plot3(L(2), 0, 0, 'rx', 'linewidth', 1.5);
% xlabel('x(n.d.)')
% ylabel('y(n.d.)')
% zlabel('z(n.d.)')
% 
% figure(70+k)
%         hold on
%         h1 = stairs(t_s(1:end)*t_star/60/60/24, [u_feasible(1,:), u_feasible(1,end)]/(1/1000/(m_0_d*l_star/t_star^2)), 'k', 'linewidth', 1.5);
%         h2 = stairs(t_s(1:end)*t_star/60/60/24, [u_feasible(1,:).*cos(u_feasible(2,:)).*cos(u_feasible(3,:)), u_feasible(1,end).*cos(u_feasible(2,end)).*cos(u_feasible(3,end))]/(1/1000/(m_0_d*l_star/t_star^2)), 'r', 'linewidth', 1.5);
%         h3 = stairs(t_s(1:end)*t_star/60/60/24, [u_feasible(1,:).*sin(u_feasible(2,:)).*cos(u_feasible(3,:)), u_feasible(1,end).*sin(u_feasible(2,end)).*cos(u_feasible(3,end))]/(1/1000/(m_0_d*l_star/t_star^2)), 'g', 'linewidth', 1.5);
%         h4 = stairs(t_s(1:end)*t_star/60/60/24, [u_feasible(1,:).*sin(u_feasible(3,:)), u_feasible(1,end).*sin(u_feasible(3,end))]/(1/1000/(m_0_d*l_star/t_star^2)), 'b', 'linewidth', 1.5);
%         xlabel('TOF(days)')
%         ylabel('Thrust(N)')
%         legend([h1, h2, h3, h4], 'T', 'Tx', 'Ty', 'Tz')
        
%% checking if any violations occur in the boundaries and linear inequalities
% 생각해보니까 어차피?



















