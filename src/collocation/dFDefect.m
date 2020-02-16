function [dFx dFu] = dFDefect(mu, Y, u, t, dt, Collocation, ...
	ispND, g0ND)

dFx = nan(21, 28);
dFu = nan(21, 3);

h = 1e-7;

Yvec = reshape(Y, [28, 1]);

for i = 1:28
	hvec = zeros(28, 1);
	hvec(i) = h;
	YvecP = Yvec + hvec;
	YvecM = Yvec - hvec;
	YP = reshape(YvecP, [7, 4]);
	YM = reshape(YvecM, [7, 4]);
	FxP = defectConstraint(mu, YP, u, t, dt, Collocation, ispND, g0ND);
	FxM = defectConstraint(mu, YM, u, t, dt, Collocation, ispND, g0ND);
	dFx(:,i) = (FxP-FxM)/2/h;
end


for i = 1:3
	
	uhvec = zeros(3, 1);
	uhvec(i) = h;
	uP = u + uhvec;
	uM = u - uhvec;
	FuP = defectConstraint(mu, Y, uP, t, dt, Collocation, ispND, g0ND);
	FuM = defectConstraint(mu, Y, uM, t, dt, Collocation, ispND, g0ND);
	dFu(:, i) = (FuP-FuM)/2/h;
end

% dFx
end