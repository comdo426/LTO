%{
Last modified: 10/19/2019
Created by: Beom Park
Name: CR3BP_LT_CSI_inertial
Description
t: now it's time-dependent
u1: control input, now constant in the inertial frame
Isp: n.d.
g0: n.d.
t0: initial time where the inertial frame and the rotational frame
coincide(usually it's assumed to be 0, but it can be changed)
%}

function dYdt = CR3BP_LT_CSI_inertial(t, Y, u1, mu, Isp, g0, t0)


    x = Y(1, :);
    y = Y(2, :);
    z = Y(3,:);
    xdot = Y(4, :);
    ydot = Y(5, :);
    zdot = Y(6, :);
    m = Y(7, :);
    T = u1(1, :);
    alpha = u1(2, :);
    alpha = alpha - (t-t0);
    beta = u1(3, :);
    
    d = sqrt((x+mu).^2 + y.^2 + z.^2);
    r = sqrt((x -1 + mu).^2 +y.^2 + z.^2);
    
%     C = [cos(t-t0), sin(t-t0), 0;
%         -sin(t-t0), cos(t-t0), 0;
%         0, 0, 1];
    
%     uhat_inert = [cos(alpha)*cos(beta); sin(alpha)*cos(beta); sin(beta)];
%     uhat_rot = C*uhat_inert;
%     uhat_rot = uhat_rot/norm(uhat_rot);
    
%     A1 = 2.*ydot + x - (1-mu)./d.^3.*(x+mu) - mu./r.^3.*(x-1+mu) + T./m.*uhat_rot(1);
%     A2 = -2.*xdot + y - (1-mu)./d.^3.*y - mu./r.^3.*y + T./m.*sin(alpha).*uhat_rot(2);
%     A3 = - (1-mu)./d.^3.*z - mu./r.^3.*z + T./m.*uhat_rot(3);

    A1 = 2.*ydot + x - (1-mu)./d.^3.*(x+mu) - mu./r.^3.*(x-1+mu) + T./m.*cos(alpha).*cos(beta);
    A2 = -2.*xdot + y - (1-mu)./d.^3.*y - mu./r.^3.*y + T./m.*sin(alpha).*cos(beta);
    A3 = - (1-mu)./d.^3.*z - mu./r.^3.*z + T./m.*sin(beta);
        
    Vel = [xdot; ydot; zdot];
    Acc = [A1; A2; A3];
    mdot = -T/Isp/g0;
    
    dYdt = [Vel; Acc; mdot];


end

