function Spacecraft = setSpacecraft(System, SC)
%SETSPACECRAFT - nondimensionlize and stores spacecraft thrust information

nPhase = length(System);
g0D = 9.80665/1000; %km/s^2
Spacecraft = cell(nPhase, 1);

for iPhase = 1:nPhase
	tstar = System{iPhase}.parameter.tstar;
	lstar = System{iPhase}.parameter.lstar;
	Spacecraft{iPhase}.thrustMaxND = SC.thrustMaxD/1000/(SC.m0D*lstar/tstar^2);
	Spacecraft{iPhase}.ispND = SC.ispD/tstar;
	Spacecraft{iPhase}.g0ND = g0D/(lstar/tstar^2);
	Spacecraft{iPhase}.thrustMaxD = SC.thrustMaxD;
	Spacecraft{iPhase}.ispD = SC.ispD;
end

end