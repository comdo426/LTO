function F = defectConstraint2BP(mu, Y, u, ~, dt, Collocation, isp, g0)

% tauRatio = Collocation.tauRatio;
phi = Collocation.phi;
phiPrime = Collocation.phiPrime;

T = u(1)*ones(1, 4);
alpha = u(2)*ones(1, 4);
beta = u(3)*ones(1, 4);

Ydot = getDerivCSI2BP(mu, Y, T, alpha, beta, isp, g0);
B = [Y, dt/2*Ydot];

xCenter = B*phi';
TCenter = u(1)*ones(1, 3);
alphaCenter = u(2)*ones(1, 3);
betaCenter = u(3)*ones(1, 3);

% Actual dynamics at even number nodes
xCenterDer = getDerivCSI2BP(mu, xCenter, TCenter, alphaCenter, betaCenter, isp, g0);
% Interpolated dynamics at even number nodes
xCenterDot = B*phiPrime';

FMatrix = xCenterDot - dt/2*xCenterDer;
F = reshape(FMatrix, [21, 1]);

end