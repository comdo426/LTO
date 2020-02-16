function dYdt = conicLT(~, Y, u1, Isp, g0)


x = Y(1);
    y = Y(2);
    z = Y(3);
    xdot = Y(4);
    ydot = Y(5);
    zdot = Y(6);
    m = Y(7);
    T = u1(1);
    alpha = u1(2);
    beta = u1(3);
	 
	 r = sqrt(x.^2 + y.^2 + z.^2);
	 
	 A1 = -x/r^3 + T./m.*cos(alpha).*cos(beta);
	 A2 = -y/r^3 + T./m.*sin(alpha).*cos(beta);
	 A3 = -z/r^3 + T./m.*sin(beta);
	 
	 mdot = -T/Isp/g0;
	 
	 dYdt = [xdot; ydot; zdot; A1; A2; A3; mdot];


end