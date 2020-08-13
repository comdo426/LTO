mu = System{1}.parameter.mu;
lstar = System{1}.parameter.lstar;
tstar = System{1}.parameter.tstar;

opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);

Y0 = InitialGuess{1}.initialConstraint(1:6);
period = InitialGuess{1}.period;
[T, Y] = ode113(@(t,y) CR3BP(t,y,mu,1,0), [0, period], Y0, opts);

figure(1)

plot3(Y(:,1), Y(:,2), Y(:,3), 'k');

JD = jday(2025,1,1,0,0,0);


mu1 = SystemPhase.parameter.mu1;
mu2 = SystemPhase.parameter.mu2;
EM.P1 = 'EARTH';
EM.P2 = 'MOON';

EM.mu = mu;
EM.lstar = lstar;
EM.tstar = tstar;
EM.mu1 = mu1;
EM.mu2 = mu2;
EM.frame = 'J2000';
EM.centralBody = 'EARTH';