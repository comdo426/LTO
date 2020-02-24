% testMergePhase
% Merges the converged solution into a single phase.

clear;

load('L2Halo2L1Halo.mat')
State = Result;
clearvars -except System State Spacecraft

% calls each phase to merge the states

[state1Mat, stateSegment1Mat, control1Mat] = getStateControlMat(State{1});
[state2Mat, stateSegment2Mat, control2Mat] = getStateControlMat(State{2});

stateSegmentMerged = [stateSegment1Mat; stateSegment2Mat(2:end, :)];

nSegment = State{1}.nSegment + State{2}.nSegment;
length(stateSegmentMerged); % 118
nSegment; % 117

timeMerged = [State{1}.timeSegment; State{1}.timeSegment(end) + State{2}.timeSegment(2:end)]; % 118

tJump = State{1}.timeSegment(end);

control2Mat(:,2) = control2Mat(:,2) + tJump;

controlMerged = [control1Mat; control2Mat];

% tstar = System{1}.parameter.tstar;
% tday = tstar/60/60/24*timeMerged;
% 
% [thrustMaxND, ~, ~, thrustMaxD] = getSpacecraftInfo(Spacecraft{1});
% 
% Ttotal = controlMerged(:,1)*thrustMaxD/thrustMaxND;
% Tx = Ttotal.*cos(controlMerged(:,2)).*cos(controlMerged(:,3));
% Ty = Ttotal.*sin(controlMerged(:,2)).*cos(controlMerged(:,3));
% Tz = Ttotal.*sin(controlMerged(:,3));
% 
% figure(3)
% hold on
% h1 = stairs(tday, [Ttotal; Ttotal(end)], 'k', 'linewidth', 1.5);
% h2 = stairs(tday, [Tx; Tx(end)], 'r', 'linewidth', 1.5);
% h3 = stairs(tday, [Ty; Ty(end)], 'g', 'linewidth', 1.5);
% h4 = stairs(tday, [Tz; Tz(end)], 'b', 'linewidth', 1.5);

% initial angle of Moon in Earth Eclip J2000 frame

load('TEST')

secPastJ2000 = (initialJDFix-cspice_j2000)*60*60*24;

moonState = cspice_spkezr('MOON', secPastJ2000, 'ECLIPJ2000', 'NONE', 'EARTH');

% figure(11)
% hold on
% axis equal
% plot3(0, 0, 0, 'bo', 'linewidth', 3.0)
% plot3(moonState(1), moonState(2), moonState(3), 'ko', 'linewidth', 2.0);

moonPos = moonState(1:3);
xUnit = [1, 0, 0];
theta0 = acos(dot(moonPos, xUnit)/norm(moonPos)/norm(xUnit));

controlMerged(:,2) = controlMerged(:,2) + theta0;

Ttotal = controlMerged(:,1)*thrustMaxD/thrustMaxND;
Tx = Ttotal.*cos(controlMerged(:,2)).*cos(controlMerged(:,3));
Ty = Ttotal.*sin(controlMerged(:,2)).*cos(controlMerged(:,3));
Tz = Ttotal.*sin(controlMerged(:,3));

figure(4)
hold on
h1 = stairs(tday, [Ttotal; Ttotal(end)], 'k', 'linewidth', 1.5);
h2 = stairs(tday, [Tx; Tx(end)], 'r', 'linewidth', 1.5);
h3 = stairs(tday, [Ty; Ty(end)], 'g', 'linewidth', 1.5);
h4 = stairs(tday, [Tz; Tz(end)], 'b', 'linewidth', 1.5);


