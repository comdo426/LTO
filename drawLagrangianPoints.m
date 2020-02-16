function lpPlot = drawLagrangianPoints(mu, pointArray)
L = Lagrange(mu);
if find(pointArray == 1) % L1
   lpPlot = plot3(L(1), 0, 0, 'rx', 'linewidth', 1.5);
end

if find(pointArray == 2) % L2
   plot3(L(2), 0, 0, 'rx', 'linewidth', 1.5);
end

if find(pointArray == 3) % L3
   plot3(L(3), 0, 0, 'rx', 'linewidth', 1.5);
end

if find(pointArray == 4) % L4
   plot3(-mu+1/2, sqrt(3)/2, 0, 'rx', 'linewidth', 1.5);
end

if find(pointArray == 5) % L5
   plot3(-mu+1/2, -sqrt(3)/2, 0, 'rx', 'linewidth', 1.5);
end
end