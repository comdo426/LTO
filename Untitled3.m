% CR3BP + 2BP 결과 분석

load('Result')

JDf = 2.460555500000000e+06 + 6.239528588041052*3.751892968837575e+05/60/60/24;

opts = odeset('AbsTol', 1e-16, 'RelTol', 3e-14);

tf = Result{2}.timeSegment(2);
x0Inert = Result{2}.state(1:7);
u1 = Result{2}.control(1:3);

isp = Spacecraft{2}.ispND;
g0 = Spacecraft{2}.g0ND;

[T, stateInert] = ode45(@(t,y) conicLT(t,y,u1, isp, g0), [0, tf], ;

figure(40)
plot3(stateInert(:,1), stateInert(:,2), stateInert(:,3))