function dF = dFInterPhaseCont(phaseFin, System1, System2, JD)

dF = zeros(7,7);
h = 1e-7;

for ii = 1:7
	
	hvec = zeros(7, 1);
	hvec(ii) = h;
	YvecP = phaseFin + hvec;
	YvecM = phaseFin - hvec;
	phaseFinP = interPhaseCont(YvecP, System1, System2, JD);
	phaseFinM = interPhaseCont(YvecM, System1, System2, JD);
	
	dF(:,ii) = (-phaseFinP + phaseFinM)/2/h;
end

end