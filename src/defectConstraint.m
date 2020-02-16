function F = defectConstraint(mu, Y, u, t, dt, Collocation, isp, g0)

tauRatio = Collocation.tauRatio;
phi = Collocation.phi;
phiPrime = Collocation.phiPrime;

T = u(1)*ones(1, 4);
angleRotated = dt*[tauRatio(1), tauRatio(3), tauRatio(5), tauRatio(7)] + ...
	t*ones(1,4);
alpha = u(2)*ones(1, 4) - angleRotated;
beta = u(3)*ones(1, 4);

Ydot = getDerivCSI(mu, Y, T, alpha, beta, isp, g0);
B = [Y, dt/2*Ydot];

xCenter = B*phi';
TCenter = u(1)*ones(1, 3);
angleCenterRotated = dt*[tauRatio(2), tauRatio(4), tauRatio(6)] + t*ones(1,3);
alphaCenter = u(2)*ones(1, 3) - angleCenterRotated;
betaCenter = u(3)*ones(1, 3);

% Actual dynamics at even number nodes
xCenterDer = getDerivCSI(mu, xCenter, TCenter, alphaCenter, betaCenter, isp, g0);
% Interpolated dynamics at even number nodes
xCenterDot = B*phiPrime';

FMatrix = xCenterDot - dt/2*xCenterDer;
F = reshape(FMatrix, [21, 1]);

end