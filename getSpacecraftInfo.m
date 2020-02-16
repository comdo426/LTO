function [thrustMaxND, ispND, g0ND, thrustMaxD, ispD] = ...
	getSpacecraftInfo(Spacecraft)

thrustMaxND = Spacecraft.thrustMaxND;
ispND = Spacecraft.ispND;
g0ND = Spacecraft.g0ND;

if nargout > 3
	thrustMaxD = Spacecraft.thrustMaxD;
	ispD = Spacecraft.ispD;
end

end