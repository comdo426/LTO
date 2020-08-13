.% Currently only applies to CR3BP - 2BP link

function phaseFinNew = interPhaseCont(phaseFin, System1, System2, JD)


EM.P1 = 'EARTH';
EM.P2 = 'MOON';
EM.mu1 = System1.parameter.mu1;
EM.mu2 = System1.parameter.mu2;
EM.frame = 'J2000';
EM.centralBody = 'EARTH';

mass = phaseFin(7);

phaseFinInert = rot2inert(phaseFin(1:6)', 0, JD, EM);

secPastJ2000 = (JD - cspice_j2000)*24*3600;
EarthVector = cspice_spkezr('EARTH', secPastJ2000, 'J2000', 'NONE', 'SUN');

phaseFinNew = EarthVector + phaseFinInert';

lstar = System2.parameter.lstar;
tstar = System2.parameter.tstar;

phaseFinNew(1:3) = phaseFinNew(1:3)/lstar;
phaseFinNew(4:6) = phaseFinNew(4:6)/lstar*tstar;

phaseFinNew = [phaseFinNew; mass];

end