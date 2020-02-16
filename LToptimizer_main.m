function Result = LToptimizer_main(SystemParameter, InitialGuess, ManifoldsParameter, option, optoption)

%% Initialization

opts = odeset('reltol', 2.3e-14, 'abstol', 1e-19);

% Get the options

EngineType = LToptget(option, 'EngineType');
FeasibilitySolver = LToptget(option, 'FeasibilitySolver');
Optimizer = LToptget(option, 'Optimizer');
MeshRefinement = LToptget(option, 'MeshRefinement');
Ephemeris = LToptget(option, 'Ephemeris');
NaturalDynamics = LToptget(option, 'NaturalDynamics');

% Set the system parameters

l_star = SystemParameter.l_star;
t_star = SystemParameter.t_star;
mu = SystemParameter.mu;
L = Lagrange(mu);
m_0_d = SystemParameter.m_0_d;
Isp_d = SystemParameter.Isp_d;
g0 =  9.80665/1000; %km/s^2

Isp_nd = Isp_d/t_star;    % Isp nondimensionalization
m_0_nd = 1; % m_0 nondimensionalization
g0_nd = g0/(l_star/t_star^2); % g0 nondimensionalization

if EngineType == 'CSI'
	T_max_d = SystemParameter.T_max_d;
	T_max = T_max_d/1000/(m_0_d*l_star/t_star^2);
else if EngineType == 'VSI'
		P_max_d = SystemParameter.P_max_d;
		P_max = P_max_d*t_star^3/m_0_d/(l_star*1000)^2;
	end
end

% Set Initial Guess parameters

x_DC_initial = InitialGuess.x;
x0_initial = InitialGuess.x0_initial;
xf_initial = InitialGuess.xf_initial;
s = InitialGuess.s;
n = 3*s;
t_s = InitialGuess.t_s;

% Set Manifolds paramters if you are using tau-alpha method

if strcmp(NaturalDynamics, 'none')
	coast = 0;
else if strcmp(NaturalDynamics, 'tau')
		coast = 2;
	else if strcmp(NaturalDynamics, 'tau-alpha')
			coast = 4;
			d = ManifoldsParameter.d;
			vu_0 = ManifoldsParameter.vu;
			vs_0 = ManifoldsParameter.vs;
			unstable_direction = ManifoldsParameter.ud;
			stable_direction = ManifoldsParameter.sd;
			P0 = ManifoldsParameter.P0;
			Pf = ManifoldsParameter.Pf;
		end
	end
end

%% Plot the initial guess

x_initial = reshape(x_DC_initial(1:21*s+7), [7, n+1]);
for i = 1:s+1
	x_initial_segment(:,i) = x_DC_initial(21*(i-1) +1: 21*(i-1)+7);
end
figure(11)
hold on
p_initial_guess = plot3(x_initial(1,:), x_initial(2,:), x_initial(3,:), 'b', 'linewidth', 1.5);
plot3(x_initial_segment(1,:), x_initial_segment(2,:), x_initial_segment(3,:), 'k.', 'linewidth', 1.0);
plot3(L(1), 0, 0, 'rx', 'linewidth', 1.5);
plot3(L(2), 0, 0, 'rx', 'linewidth', 1.5);
xlabel('x(n.d.)')
ylabel('y(n.d.)')
zlabel('z(n.d.)')
axis equal

if strcmp(NaturalDynamics, 'tau-alpha')
	[T0 Y0] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 1), [0 x_DC_initial(21*s+7+3*s+1)], [x0_initial; reshape(eye(6), [36, 1])], opts);
	[Tf Yf] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 1), [0 x_DC_initial(21*s+7+3*s+4)], [xf_initial; reshape(eye(6), [36, 1])], opts);
	xu = Y0(end, 1:6)';
	xs = Yf(end, 1:6)';
	
	phi_0 = reshape(Y0(end, 7:42), [6, 6]);
	phi_f = reshape(Yf(end, 7:42), [6, 6]);
	
	vu = phi_0*vu_0;
	vu = vu/norm(vu);
	vs = phi_f*vs_0;
	vs = vs/norm(vs);
	
	xu_step = xu + unstable_direction*d*vu;
	xs_step = xs + stable_direction*d*vs;
	
	[Tu Yu] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), [0 x_DC_initial(21*s+7+3*s+2)], xu_step, opts);
	[Ts Ys] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), [0 x_DC_initial(21*s+7+3*s+3)], xs_step, opts);
	
	[Tf Yf] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), [0 Pf-x_DC_initial(21*s+7+3*s+4)], xs, opts);
	
	figure(11)
	hold on
	p_initial1 = plot3(Y0(:, 1), Y0(:,2), Y0(:,3), '--k', 'linewidth', 1.5)
	p_initial2 = plot3(Yf(:,1),Yf(:,2),Yf(:,3), '--k', 'linewidth', 1.5)
	p_initial3 = plot3(Yu(:, 1), Yu(:,2), Yu(:,3), 'k', 'linewidth', 1.5)
	p_initial4 = plot3(Ys(:, 1), Ys(:,2), Ys(:,3), 'k', 'linewidth', 1.5)
	DrawMoon(11, mu, l_star, 1);
	
	legend([p_initial_guess, p_initial1, p_initial3], 'Initial guess(collocation)', 'periodic orbits', 'manifolds')
end
% DrawEarth(101, mu, l_star, 1);


%% Initialization for the solvers

[Phi, Phi_prime, Phi_mesh_add] = LGL_7th_coefficient;
tau = [-1; -sqrt(495 +66*sqrt(15))/33; -sqrt(495 -66*sqrt(15))/33; 0; sqrt(495 -66*sqrt(15))/33; sqrt(495 +66*sqrt(15))/33; 1];

k = 1;

OUTPUT.done_switch = 0;
OUTPUT.remove_switch = 1;

while OUTPUT.done_switch ~= 1
	
	if FeasibilitySolver == 1
		if k < 1 % if you are solving it for the first time (새로 푸는 거면 이거 k == 1로 바꿔야함)
			x_DC_initial = [x_DC_initial; 0.0001*ones(s,1)];
		else
			for i = 1:s
				x_DC_initial = [x_DC_initial; asin(sqrt(x_DC_initial(21*s+7+3*(i-1)+1)/T_max))];
			end
		end
		%
		%     x_DC_initial_feasibility = x_DC_initial;
		%
		%     x_DC_converged_feasibility = Feasibility_solver_CSI_orbit_stack(mu, x_DC_initial_feasibility,  t_s, s, n, Phi, Phi_prime, x0_initial, xf_initial, T_max, Isp_nd, g0_nd);
		
		beq =  [1];
		A = [];
		b = [];
		if EngineType == 'CSI'
			Aeq = zeros(1, 7*(n+1) + 3*s + coast); % Size of the equality constraint matrix
			Aeq(1, 7) = 1;
			lb = -inf*ones(7*(n+1) + 3*s + coast, 1);
			ub = inf*ones(7*(n+1) + 3*s + coast, 1);
			for i = 1:s
				lb(7*(n+1) + 3*(i-1) + 1) = 0; % little bit of deviation for IP method
				ub(7*(n+1) + 3*(i-1) + 1) = T_max; % little bit of deviation for IP method
			end
		else
			Aeq = zeros(1, 7*(n+1) + 4*s + coast);
			Aeq(1, 7) = 1;
			lb = -inf*ones(7*(n+1) + 4*s + coast, 1);
			ub = inf*ones(7*(n+1) + 4*s + coast, 1);
			for i = 1:s
				lb(7*(n+1) + 4*(i-1) +1) = 0-eps; % little bit of deviation for IP method
				ub(7*(n+1) + 4*(i-1) +1) = P_max + eps(P_max); % little bit of deviation for IP method
			end
		end
		%     feasible = @(x) 1;
		if k == 1
			opts2 = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display','iter', 'MaxFunctionEvaluations', 30000000, 'StepTolerance', 1e-15, 'CheckGradient', true, 'SpecifyConstraintGradient', true, 'SpecifyObjectiveGradient', true, 'FiniteDifferenceType', 'central', 'ConstraintTolerance', 1e-13, 'OptimalityTolerance', 1, 'maxiteration',3000);
		else
			opts2 = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display','iter', 'MaxFunctionEvaluations', 30000000, 'StepTolerance', 1e-15, 'CheckGradient', false, 'SpecifyConstraintGradient', true, 'SpecifyObjectiveGradient', true, 'FiniteDifferenceType', 'central', 'ConstraintTolerance', 1e-13, 'OptimalityTolerance', 1, 'maxiteration',3000);
			opts3 = optimoptions('fsolve', 'Display','iter', 'MaxFunctionEvaluations', 30000000, 'StepTolerance', 1e-15, 'CheckGradient', false, 'SpecifyObjectiveGradient', true, 'OptimalityTolerance', 1e-10, 'maxiteration',3000);
		end
		if strcmp(NaturalDynamics, 'none')
			x_DC_converged_feasibility = Feasibility_solver_CSI_orbit_stack_inertial(mu, x_DC_initial, t_s, s, n, Phi, Phi_prime, x0_initial, xf_initial, T_max, Isp_nd, g0_nd, 0);
			%         FeasibilitySolver = 0;,
			
			%         if k == 1
			%             x_DC_converged_feasibility = fmincon(@(x) feasible(x, s) , x_DC_initial, A, b, Aeq, beq, lb, ub, @(x) LGL_7th_with_control_CSI_orbit_stack(mu, x, t_s, s, n, Phi, Phi_prime, x0_initial, xf_initial, Isp_nd, g0_nd), opts2);
			%         else
			% % %                 save('XDC', 'x_DC_initial')
			% %                 size(x_DC_initial)
			% %                 for i = 1:s
			% %                     if x_DC_initial(21*s+7+3*(i-1)+1) < 0
			% %                         x_DC_initial(21*s+7+3*(i-1)+1) = eps(1);
			% %                     end
			% %                     x_DC_initial = [x_DC_initial; asin(sqrt(x_DC_initial(21*s+7+3*(i-1)+1)/T_max))];
			% %                 end
			% %                 size(x_DC_initial)
			% % %                 imag(x_DC_initial)
			% %             x_DC_converged_feasibility =  fsolve(@(x) Feasibility_solver_CSI_orbit_stack_fsolve(mu, x, t_s, s, n, Phi, Phi_prime, x0_initial, xf_initial, T_max, Isp_nd, g0_nd), x_DC_initial, opts3);
			% %             x_DC_converged_feasibility = x_DC_converged_feasibility(1:21*s +7 + 3*s);
			%             x_DC_converged_feasibility = fmincon(@(x) feasible(x, s) , x_DC_initial, A, b, Aeq, beq, lb, ub, @(x) LGL_7th_with_control_CSI_orbit_stack(mu, x, t_s, s, n, Phi, Phi_prime, x0_initial, xf_initial, Isp_nd, g0_nd), opts2);
			%         end
			
		else if strcmp(NaturalDynamics, 'tau')
				
			else if strcmp(NaturalDynamics, 'tau-alpha')
					x_DC_converged_feasibility = fmincon(@(x) feasible(x, s), x_DC_initial, A, b, Aeq, beq, lb, ub, @(x) LGL_7th_with_control_CSI_manifolds(mu, x, t_s, s, n, Phi, Phi_prime, x0_initial, xf_initial, d, unstable_direction, stable_direction, vu_0, vs_0, Isp_nd, g0_nd), opts2);
				end
			end
		end
		
		
		
		if k > 0
			%     if k == 1
			x_feasible = reshape(x_DC_converged_feasibility(1:21*s +7), [7, n+1]);
			for i = 1:s+1
				x_feasible_segment(:,i) = x_DC_converged_feasibility(21*(i-1) +1: 21*(i-1)+7);
			end
			u_feasible = reshape(x_DC_converged_feasibility(21*s + 7 +1:21*s+7+3*s), [3, s]);
			figure(50+k)
			hold on
			for i = 1:s
				if u_feasible(1,i) < 0.001
					color = [0 0 1];
				else
					color = [1 0 0];
				end
				p_feasible = plot3(x_feasible(1, 3*(i-1)+1:3*i+1), x_feasible(2, 3*(i-1)+1:3*i+1), x_feasible(3, 3*(i-1)+1:3*i+1), 'Color', color, 'linewidth', 1.5);
			end
			p4 = plot3(x_feasible_segment(1, :), x_feasible_segment(2,:), x_feasible_segment(3,:), 'k.', 'linewidth', 1.0);
			axis equal
			plot3(L(1), 0, 0, 'rx', 'linewidth', 1.5);
			plot3(L(2), 0, 0, 'rx', 'linewidth', 1.5);
			xlabel('x(n.d.)')
			ylabel('y(n.d.)')
			zlabel('z(n.d.)')
			
			if strcmp(NaturalDynamics, 'tau-alpha')
				[T0 Y0] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 1), [0 x_DC_converged_feasibility(21*s+7+3*s+1)], [x0_initial; reshape(eye(6), [36, 1])], opts);
				[Tf Yf] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 1), [0 x_DC_converged_feasibility(21*s+7+3*s+4)], [xf_initial; reshape(eye(6), [36, 1])], opts);
				xu = Y0(end, 1:6)';
				xs = Yf(end, 1:6)';
				
				phi_0 = reshape(Y0(end, 7:42), [6, 6]);
				phi_f = reshape(Yf(end, 7:42), [6, 6]);
				
				vu = phi_0*vu_0;
				vu = vu/norm(vu);
				vs = phi_f*vs_0;
				vs = vs/norm(vs);
				
				xu_step = xu + unstable_direction*d*vu;
				xs_step = xs + stable_direction*d*vs;
				
				[Tu Yu] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), [0 x_DC_converged_feasibility(21*s+7+3*s+2)], xu_step, opts);
				[Ts Ys] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), [0 x_DC_converged_feasibility(21*s+7+3*s+3)], xs_step, opts);
				
				[Tf Yf] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), [0 Pf-x_DC_converged_feasibility(21*s+7+3*s+4)], xs, opts);
				
				figure(50+k)
				hold on
				p_feasible1 = plot3(Y0(:, 1), Y0(:,2), Y0(:,3), '--k', 'linewidth', 1.5)
				p_feasible2 = plot3(Yf(:,1),Yf(:,2),Yf(:,3), '--k', 'linewidth', 1.5)
				p_feasible3 = plot3(Yu(:, 1), Yu(:,2), Yu(:,3), 'k', 'linewidth', 1.5)
				p_feasible4 = plot3(Ys(:, 1), Ys(:,2), Ys(:,3), 'k', 'linewidth', 1.5)
				DrawMoon(50+k, mu, l_star, 1);
				
				legend([p_feasible, p_feasible1, p_feasible3], 'Feasible transfer', 'periodic orbits', 'manifolds')
			end
			
			
			figure(70+k)
			hold on
			h1 = stairs(t_s(1:end)*t_star/60/60/24, [u_feasible(1,:), u_feasible(1,end)]/(1/1000/(m_0_d*l_star/t_star^2)), 'k', 'linewidth', 1.5);
			h2 = stairs(t_s(1:end)*t_star/60/60/24, [u_feasible(1,:).*cos(u_feasible(2,:)).*cos(u_feasible(3,:)), u_feasible(1,end).*cos(u_feasible(2,end)).*cos(u_feasible(3,end))]/(1/1000/(m_0_d*l_star/t_star^2)), 'r', 'linewidth', 1.5);
			h3 = stairs(t_s(1:end)*t_star/60/60/24, [u_feasible(1,:).*sin(u_feasible(2,:)).*cos(u_feasible(3,:)), u_feasible(1,end).*sin(u_feasible(2,end)).*cos(u_feasible(3,end))]/(1/1000/(m_0_d*l_star/t_star^2)), 'g', 'linewidth', 1.5);
			h4 = stairs(t_s(1:end)*t_star/60/60/24, [u_feasible(1,:).*sin(u_feasible(3,:)), u_feasible(1,end).*sin(u_feasible(3,end))]/(1/1000/(m_0_d*l_star/t_star^2)), 'b', 'linewidth', 1.5);
			xlabel('TOF(days)')
			ylabel('Thrust(N)')
			legend([h1, h2, h3, h4], 'T', 'Tx', 'Ty', 'Tz')
			Result.m_feasible = x_feasible(end);
			Result.x_feasible = x_DC_converged_feasibility;
			Result.t_s_feasible = t_s;
			Result.s_feasible = s;
			
			prompt_file = 'Feasibility solved, file save?'
			save_flag = input(prompt_file);
			if save_flag == 0
			else
				save(save_flag, 'Result')
			end
			prompt = 'Feasibility Solved, Continue?';
			Continue_flag = input(prompt);
			if Continue_flag == 1
			else
				break;
			end
		end
	else
		x_DC_converged_feasibility = x_DC_initial;
	end
	
	final_mass = @(x) -x(7*n + 7);
	% Constraining initial mass to 1 (m0 = 1)
	beq =  [1];
	A = [];
	b = [];
	if EngineType == 'CSI'
		Aeq = zeros(1, 7*(n+1) + 3*s + coast); % Size of the equality constraint matrix
		Aeq(1, 7) = 1;
		lb = -inf*ones(7*(n+1) + 3*s + coast, 1);
		ub = inf*ones(7*(n+1) + 3*s + coast, 1);
		for i = 1:s
			if x_DC_converged_feasibility(21*s + 7 + 3*(i-1) + 1) < eps
				x_DC_converged_feasibility(21*s + 7 + 3*(i-1) + 1) = 0;
			else if T_max - x_DC_converged_feasibility(21*s + 7 + 3*(i-1) + 1) < eps(T_max)
					x_DC_converged_feasibility(21*s + 7 + 3*(i-1) + 1) = T_max - eps(T_max);
				end
			end
			lb(7*(n+1) + 3*(i-1) + 1) = 0; % little bit of deviation for IP method
			ub(7*(n+1) + 3*(i-1)  +1) = T_max; % little bit of deviation for IP method
		end
	else
		Aeq = zeros(1, 7*(n+1) + 4*s + coast);
		Aeq(1, 7) = 1;
		lb = -inf*ones(7*(n+1) + 4*s + coast, 1);
		ub = inf*ones(7*(n+1) + 4*s + coast, 1);
		for i = 1:s
			lb(7*(n+1) + 4*(i-1) +1) = 0-eps; % little bit of deviation for IP method
			ub(7*(n+1) + 4*(i-1) +1) = P_max + eps(P_max); % little bit of deviation for IP method
		end
	end
	
	if Optimizer == 1
		t_s_old = t_s;
		tic
		if strcmp(NaturalDynamics, 'none')
			if k>1
				optoption = optimoptions(optoption, 'InitBarrierParam', 1e-3);
			end
			[x_DC_converged ,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x) final_mass_fcn(x, s), x_DC_converged_feasibility, A, b, Aeq, beq, lb, ub, @(x) LGL_7th_with_control_CSI_orbit_stack_inertial(mu, x, t_s, s, n, Phi, Phi_prime, x0_initial, xf_initial, Isp_nd, g0_nd, 0), optoption);
			save('output_test', 'output');
		else if strcmp(NaturalDynamics, 'tau')
				
			else if strcmp(NaturalDynamics, 'tau-alpha')
					[x_DC_converged ,fval,exitflag,output,lambda,grad,hessian] = fmincon(final_mass_fcn, x_DC_converged_feasibility, A, b, Aeq, beq, lb, ub, @(x) LGL_7th_with_control_CSI_manifolds(mu, x, t_s, s, n, Phi, Phi_prime, x0_initial, xf_initial, d, unstable_direction, stable_direction, vu_0, vs_0, Isp_nd, g0_nd), optoption);
				end
			end
		end
		toc
		%     if k == 1
		x_optimized = reshape(x_DC_converged(1:21*s +7), [7, n+1]);
		for i = 1:s+1
			x_optimized_segment(:,i) = x_DC_converged(21*(i-1) +1: 21*(i-1)+7);
		end
		u_optimized = reshape(x_DC_converged(21*s + 7 +1:21*s+7+3*s), [3, s]);
		figure(100+k)
		hold on
		for i = 1:s
			if u_optimized(1,i) < 0.0001
				color = [0 0 1];
			else
				color = [1 0 0];
			end
			p_optimized = plot3(x_optimized(1, 3*(i-1)+1:3*i+1), x_optimized(2, 3*(i-1)+1:3*i+1), x_optimized(3, 3*(i-1)+1:3*i+1), 'Color', color, 'linewidth', 1.5);
		end
		p4 = plot3(x_optimized_segment(1, :), x_optimized_segment(2,:), x_optimized_segment(3,:), 'k.', 'linewidth', 0.3);
		plot3(L(1), 0, 0, 'kx', 'linewidth', 1.5);
		plot3(L(2), 0, 0, 'kx', 'linewidth', 1.5);
		xlabel('x(n.d.)')
		ylabel('y(n.d.)')
		zlabel('z(n.d.)')
		axis equal
		
		
		if strcmp(NaturalDynamics, 'tau-alpha')
			[T0 Y0] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 1), [0 x_DC_converged(21*s+7+3*s+1)], [x0_initial; reshape(eye(6), [36, 1])], opts);
			[Tf Yf] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 1), [0 x_DC_converged(21*s+7+3*s+4)], [xf_initial; reshape(eye(6), [36, 1])], opts);
			xu = Y0(end, 1:6)';
			xs = Yf(end, 1:6)';
			
			phi_0 = reshape(Y0(end, 7:42), [6, 6]);
			phi_f = reshape(Yf(end, 7:42), [6, 6]);
			
			vu = phi_0*vu_0;
			vu = vu/norm(vu);
			vs = phi_f*vs_0;
			vs = vs/norm(vs);
			
			xu_step = xu + unstable_direction*d*vu;
			xs_step = xs + stable_direction*d*vs;
			
			[Tu Yu] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), [0 x_DC_converged(21*s+7+3*s+2)], xu_step, opts);
			[Ts Ys] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), [0 x_DC_converged(21*s+7+3*s+3)], xs_step, opts);
			
			[Tf Yf] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), [0 Pf-x_DC_converged(21*s+7+3*s+4)], xs, opts);
			
			figure(100+k)
			hold on
			p_feasible1 = plot3(Y0(:, 1), Y0(:,2), Y0(:,3), '--k', 'linewidth', 1.5)
			p_feasible2 = plot3(Yf(:,1),Yf(:,2),Yf(:,3), '--k', 'linewidth', 1.5)
			p_feasible3 = plot3(Yu(:, 1), Yu(:,2), Yu(:,3), 'k', 'linewidth', 1.5)
			p_feasible4 = plot3(Ys(:, 1), Ys(:,2), Ys(:,3), 'k', 'linewidth', 1.5)
			DrawMoon(100+k, mu, l_star, 1);
			
			%             legend([p_feasible, p_feasible1, p_feasible3], 'Feasible transfer', 'periodic orbits', 'manifolds')
		end
		
		figure(120+k)
		hold on
		h1 = stairs(t_s(1:end)*t_star/60/60/24, [u_optimized(1,:), u_optimized(1,end)]/(1/1000/(m_0_d*l_star/t_star^2)), 'k', 'linewidth', 1.5);
		h2 = stairs(t_s(1:end)*t_star/60/60/24, [u_optimized(1,:).*cos(u_optimized(2,:)).*cos(u_optimized(3,:)), u_optimized(1,end).*cos(u_optimized(2,end)).*cos(u_optimized(3,end))]/(1/1000/(m_0_d*l_star/t_star^2)), 'r', 'linewidth', 1.5);
		h3 = stairs(t_s(1:end)*t_star/60/60/24, [u_optimized(1,:).*sin(u_optimized(2,:)).*cos(u_optimized(3,:)), u_optimized(1,end).*sin(u_optimized(2,end)).*cos(u_optimized(3,end))]/(1/1000/(m_0_d*l_star/t_star^2)), 'g', 'linewidth', 1.5);
		h4 = stairs(t_s(1:end)*t_star/60/60/24, [u_optimized(1,:).*sin(u_optimized(3,:)), u_optimized(1,end).*sin(u_optimized(3,end))]/(1/1000/(m_0_d*l_star/t_star^2)), 'b', 'linewidth', 1.5);
		xlabel('TOF(days)')
		ylabel('Thrust(N)')
		legend([h1, h2, h3, h4], 'T', 'Tx', 'Ty', 'Tz')
		Result.m_optimized = x_optimized(end);
		Result.x_optimized = x_DC_converged;
		Result.t_s_feasible = t_s;
		Result.s_feasible = s;
		
		prompt_file = 'Optimizer solved, file save?'
		save_flag = input(prompt_file);
		if save_flag == 0
		else
			save(save_flag, 'Result')
		end
		prompt = 'Optimizer Solved, Continue?';
		Continue_flag = input(prompt);
		if Continue_flag == 1
		else
			break;
		end
		%     end
	else
		x_DC_converged = x_DC_converged_feasibility;
	end
	
	if MeshRefinement ~= 1
		break;
	end
	
	if OUTPUT.remove_switch == 1
		sprintf('Removing unnecessary segments...: %d', k)
	else
		sprintf('Adding necessary segments...: %d', k)
	end
	
	Result.m_mesh(k) = x_DC_converged(7*n+7);
	
	OUTPUT = CEP_mesh_refinement_CSI_inertial(mu, x_DC_converged, s, t_s, OUTPUT.remove_switch, Isp_nd, g0_nd, Phi_mesh_add);
	s = OUTPUT.s;
	n = 3*s;
	t_s = OUTPUT.time;
	x_DC_initial = OUTPUT.x;
	
	x_DC_mesh = x_DC_initial;
	%     save('before_mesh.mat', 'x_DC_converged');
	%     save('after_mesh.mat', 'x_DC_mesh');
	
	%     prompt = 'Mesh Refinement Step: Continue?';
	%         Continue_flag = input(prompt);
	%         if Continue_flag == 1
	%         else
	%             break;
	%         end
	prompt = 'Mesh Refinement visual inspection?';
	Mesh_Refinement_visual_inspection = input(prompt);
	
	if Mesh_Refinement_visual_inspection == 1
		x_optimized_refined = reshape(x_DC_mesh(1:21*s +7), [7, n+1]);
		u_optimized_refined = reshape(x_DC_mesh(21*s + 7 +1:21*s+7+3*s), [3, s]);
		figure(200 + k) % position graph
		hold on
		plot3(x_optimized_refined(1,:), x_optimized_refined(2,:), x_optimized_refined(3,:), '-r', 'linewidth', 1.5)
		plot3(x_optimized(1,:), x_optimized(2,:), x_optimized(3,:), '-k', 'linewidth', 1.5)
		plot3(x_optimized_refined(1,:), x_optimized_refined(2,:), x_optimized_refined(3,:), 'rx', 'linewidth', 1.5)
		plot3(x_optimized(1,:), x_optimized(2,:), x_optimized(3,:), 'kx', 'linewidth', 1.5)
		figure(220 + k) % velocity graph
		hold on
		plot3(x_optimized_refined(4,:), x_optimized_refined(5,:), x_optimized_refined(6,:), '-r', 'linewidth', 1.5)
		plot3(x_optimized(4,:), x_optimized(5,:), x_optimized(6,:), '-k', 'linewidth', 1.5)
		%         figure(240 + k) % mass graph
		%             plot(t_s
		figure(260 + k) % thrust graph
		hold on
		plot(t_s(1:end-1), u_optimized_refined(1, :), 'rx')
		plot(t_s_old(1:end-1), u_optimized(1,:), '-k')
		figure(280 + k) % thrust graph
		hold on
		plot(t_s(1:end-1), u_optimized_refined(2, :), 'rx')
		plot(t_s_old(1:end-1), u_optimized(2,:), '-k')
		figure(300 + k) % thrust graph
		hold on
		plot(t_s(1:end-1), u_optimized_refined(3, :), 'rx')
		plot(t_s_old(1:end-1), u_optimized(3,:), '-k')
	else
		
	end
	
	
	
	k = k+1;
	
	fprintf('Mesh done: %d\n', OUTPUT.done_switch);
	
	
	
	
	
	
	
end

%
%     if MeshRefinement == 1
%         x_mesh = reshape(x_DC_mesh(1:21*s +7), [7, n+1]);
%         for i = 1:s+1
%             x_mesh_segment(:,i) = x_DC_mesh(21*(i-1) +1: 21*(i-1)+7);
%         end
%         u_mesh = reshape(x_DC_mesh(21*s + 7 +1:21*s+7+3*s), [3, s]);
%         figure(150)
%         hold on
%         for i = 1:s
%             if u_mesh(1,i) < 0.01
%                 color = [0 0 1];
%             else
%                 color = [1 0 0];
%             end
%             p4 = plot3(x_mesh(1, 3*(i-1)+1:3*i+1), x_mesh(2, 3*(i-1)+1:3*i+1), x_mesh(3, 3*(i-1)+1:3*i+1), 'Color', color, 'linewidth', 1.5);
%         end
%         p4 = plot3(x_mesh_segment(1, :), x_mesh_segment(2,:), x_mesh_segment(3,:), 'k.', 'linewidth', 0.3);
%         plot3(L(1), 0, 0, 'kx', 'linewidth', 1.5);
%         plot3(L(2), 0, 0, 'kx', 'linewidth', 1.5);
%         xlabel('x(n.d.)')
%         ylabel('y(n.d.)')
%         zlabel('z(n.d.)')
%         axis equal
%
%         if strcmp(NaturalDynamics, 'tau-alpha')
%             [T0 Y0] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 1), [0 x_DC_mesh(21*s+7+3*s+1)], [x0_initial; reshape(eye(6), [36, 1])], opts);
%             [Tf Yf] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 1), [0 x_DC_mesh(21*s+7+3*s+4)], [xf_initial; reshape(eye(6), [36, 1])], opts);
%             xu = Y0(end, 1:6)';
%             xs = Yf(end, 1:6)';
%
%             phi_0 = reshape(Y0(end, 7:42), [6, 6]);
%             phi_f = reshape(Yf(end, 7:42), [6, 6]);
%
%             vu = phi_0*vu_0;
%             vu = vu/norm(vu);
%             vs = phi_f*vs_0;
%             vs = vs/norm(vs);
%
%             xu_step = xu + unstable_direction*d*vu;
%             xs_step = xs + stable_direction*d*vs;
%
%             [Tu Yu] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), [0 x_DC_mesh(21*s+7+3*s+2)], xu_step, opts);
%             [Ts Ys] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), [0 x_DC_mesh(21*s+7+3*s+3)], xs_step, opts);
%
%             [Tf Yf] = ode113(@(t,y) CR3BP_master(t, y, mu, 1, 0), [0 Pf-x_DC_mesh(21*s+7+3*s+4)], xs, opts);
%
%             figure(150)
%             hold on
%             p_feasible1 = plot3(Y0(:, 1), Y0(:,2), Y0(:,3), '--k', 'linewidth', 1.5)
%             p_feasible2 = plot3(Yf(:,1),Yf(:,2),Yf(:,3), '--k', 'linewidth', 1.5)
%             p_feasible3 = plot3(Yu(:, 1), Yu(:,2), Yu(:,3), 'k', 'linewidth', 1.5)
%             p_feasible4 = plot3(Ys(:, 1), Ys(:,2), Ys(:,3), 'k', 'linewidth', 1.5)
%             DrawMoon(150, mu, l_star, 1);
%
% %             legend([p_feasible, p_feasible1, p_feasible3], 'Feasible transfer', 'periodic orbits', 'manifolds')
%         end
%
% %         DrawMoon(106, mu, l_star, 10);
%         figure(171)
%         hold on
%         h1 = stairs(t_s(1:end)*t_star/60/60/24, [u_mesh(1,:), u_mesh(1,end)]/(1/1000/(m_0_d*l_star/t_star^2)), 'k', 'linewidth', 1.5);
%         h2 = stairs(t_s(1:end)*t_star/60/60/24, [u_mesh(1,:).*cos(u_mesh(2,:)).*cos(u_mesh(3,:)), u_mesh(1,end).*cos(u_mesh(2,end)).*cos(u_mesh(3,end))]/(1/1000/(m_0_d*l_star/t_star^2)), 'r', 'linewidth', 1.5);
%         h3 = stairs(t_s(1:end)*t_star/60/60/24, [u_mesh(1,:).*sin(u_mesh(2,:)).*cos(u_mesh(3,:)), u_mesh(1,end).*sin(u_mesh(2,end)).*cos(u_mesh(3,end))]/(1/1000/(m_0_d*l_star/t_star^2)), 'g', 'linewidth', 1.5);
%         h4 = stairs(t_s(1:end)*t_star/60/60/24, [u_mesh(1,:).*sin(u_mesh(3,:)), u_mesh(1,end).*sin(u_mesh(3,end))]/(1/1000/(m_0_d*l_star/t_star^2)), 'b', 'linewidth', 1.5);
%
%         xlabel('TOF(days)')
%         ylabel('Thrust(N)')
%         legend([h1, h2, h3, h4], 'T', 'Tx', 'Ty', 'Tz')
%         %         break
%         Result.OUTPUT = OUTPUT;
%         Result.x_mesh = x_DC_initial;
%     end

end