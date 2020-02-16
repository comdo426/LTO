function [c F GC dF] = LGL_7th_with_control_CSI_orbit_stack(mu, x,  t, s, n, Phi, Phi_prime, x0_initial, xf_initial,  Isp_nd, g0_nd, t0)
    opts = odeset('reltol', 2.3e-14, 'abstol', 1e-19);
    
    cnt = 0;
    tau = [-1; -sqrt(495 +66*sqrt(15))/33; -sqrt(495 -66*sqrt(15))/33; 0; sqrt(495 -66*sqrt(15))/33; sqrt(495 +66*sqrt(15))/33; 1];
    
    c = [];
        cnt = 0;
    for k = 1:s
%         if k == 1
%             x1(:, k) = [x0_initial; 1];
%         else
            x1(:, k) = x(28*(k-1)+1-7*cnt:28*(k-1)+7-7*cnt);
%         end
            x3(:, k) = x(28*(k-1)+8-7*cnt:28*(k-1)+14-7*cnt);
            x5(:, k) = x(28*(k-1)+15-7*cnt:28*(k-1)+21-7*cnt);
%         if k == s
%             x7(1:6, k) = xf_initial;
%             x7(7, k) = x(28*(k-1)+28-7*cnt);
%         else
            x7(:, k) = x(28*(k-1)+22-7*cnt:28*(k-1)+28-7*cnt);
%         end
            u1(:, k) = x(7*(n+1) + 3*(k-1)+1: 7*(n+1) + 3*k);
            u3(:, k) = x(7*(n+1) + 3*(k-1)+1: 7*(n+1) + 3*k);
            u5(:, k) = x(7*(n+1) + 3*(k-1)+1: 7*(n+1) + 3*k);
            u7(:, k) = x(7*(n+1) + 3*(k-1)+1: 7*(n+1) + 3*k);
            cnt = cnt + 1;
    end
    
    
    for k = 1:s
        Delta_t(k) = t(k+1) - t(k);
    end
    
    
    for k = 1:s
           
        xdot1(:, k) = GetVelAcc_6dim_LT_CSI_inertial(t(k), mu, x1(:,k), u1(:,k), Isp_nd, g0_nd, 0);
        xdot3(:, k) = GetVelAcc_6dim_LT_CSI_inertial(t(k) + Delta_t(k)*(tau(3)+1)/2, mu, x3(:,k), u3(:,k), Isp_nd, g0_nd, 0);    
        xdot5(:, k) = GetVelAcc_6dim_LT_CSI_inertial(t(k) + Delta_t(k)*(tau(5)+1)/2, mu, x5(:,k), u5(:,k), Isp_nd, g0_nd, 0);
        xdot7(:, k) = GetVelAcc_6dim_LT_CSI_inertial(t(k+1), mu, x7(:,k), u7(:,k), Isp_nd, g0_nd, 0);
        
        B = [x1(:, k), x3(:, k), x5(:, k), x7(:, k), Delta_t(k)/2*xdot1(:,k), Delta_t(k)/2*xdot3(:,k), Delta_t(k)/2*xdot5(:,k), Delta_t(k)/2*xdot7(:,k)];
        
        xc = B*Phi';
        x2(:,k) = xc(:, 1);
        x4(:,k) = xc(:, 2);
        x6(:,k) = xc(:, 3);
        
        u2(:,k) = u1(:,k);
        u4(:,k) = u3(:,k);
        u6(:,k) = u5(:,k);
        
        f2(:,k) = GetVelAcc_6dim_LT_CSI_inertial(t(k) + Delta_t(k)*(tau(2)+1)/2, mu, x2(:,k), u2(:,k), Isp_nd, g0_nd, 0);
        f4(:,k) = GetVelAcc_6dim_LT_CSI_inertial(t(k) + Delta_t(k)*(tau(4)+1)/2, mu, x4(:,k), u4(:,k), Isp_nd, g0_nd, 0);
        f6(:,k) = GetVelAcc_6dim_LT_CSI_inertial(t(k) + Delta_t(k)*(tau(6)+1)/2, mu, x6(:,k), u6(:,k), Isp_nd, g0_nd, 0);
        
        xcdot = B*Phi_prime';
              
        Del_2 = xcdot(:,1) - Delta_t(k)/2*f2(:,k);
        Del_4 = xcdot(:,2) - Delta_t(k)/2*f4(:,k);
        Del_6 = xcdot(:,3) - Delta_t(k)/2*f6(:,k);
        
        F(21*(k-1)+1: 21*k, 1) = [Del_2; Del_4; Del_6];
    end

    F(21*s + 1 : 21*s+6) = x0_initial - x1(1:6,1);
    F(21*s + 7 : 21*s +12) = xf_initial - x7(1:6, end);
    
    
    if nargout > 1
        
        
        GC = [];
        dudu = reshape(eye(3), [9, 1]);
        
        cnt = 0;
        
        for k = 1:s
            
            C1 = [cos(t(k)-t0), sin(t(k)-t0), 0;
                -sin(t(k)-t0), cos(t(k)-t0), 0;
                0, 0, 1];
            C2 = [cos(t(k)+Delta_t(k)/2*(tau(2)+1)-t0), sin(t(k)+Delta_t(k)/2*(tau(2)+1)-t0), 0;
                -sin(t(k)+Delta_t(k)/2*(tau(2)+1)-t0), cos(t(k)+Delta_t(k)/2*(tau(2)+1)-t0), 0;
                0, 0, 1];
            C3 = [cos(t(k)+Delta_t(k)/2*(tau(3)+1)-t0), sin(t(k)+Delta_t(k)/2*(tau(3)+1)-t0), 0;
                -sin(t(k)+Delta_t(k)/2*(tau(3)+1)-t0), cos(t(k)+Delta_t(k)/2*(tau(3)+1)-t0), 0;
                0, 0, 1];
            C4 = [cos(t(k)+Delta_t(k)/2*(tau(4)+1)-t0), sin(t(k)+Delta_t(k)/2*(tau(4)+1)-t0), 0;
                -sin(t(k)+Delta_t(k)/2*(tau(4)+1)-t0), cos(t(k)+Delta_t(k)/2*(tau(4)+1)-t0), 0;
                0, 0, 1];
            C5 = [cos(t(k)+Delta_t(k)/2*(tau(5)+1)-t0), sin(t(k)+Delta_t(k)/2*(tau(5)+1)-t0), 0;
                -sin(t(k)+Delta_t(k)/2*(tau(5)+1)-t0), cos(t(k)+Delta_t(k)/2*(tau(5)+1)-t0), 0;
                0, 0, 1];
            C6 = [cos(t(k)+Delta_t(k)/2*(tau(6)+1)-t0), sin(t(k)+Delta_t(k)/2*(tau(6)+1)-t0), 0;
                -sin(t(k)+Delta_t(k)/2*(tau(6)+1)-t0), cos(t(k)+Delta_t(k)/2*(tau(6)+1)-t0), 0;
                0, 0, 1];
            C7 = [cos(t(k+1)-t0), sin(t(k+1)-t0), 0;
                -sin(t(k+1)-t0), cos(t(k+1)-t0), 0;
                0, 0, 1];
            
            uhat_inert = [cos(u1(2,k))*cos(u1(3,k)); sin(u1(2,k))*cos(u1(3,k)); sin(u1(3,k))];
            
            [U_xx1(k) U_xy1(k) U_xz1(k) U_yy1(k) U_yz1(k) U_zz1(k)] = GetDerivative_6dim(mu, x1(:,k));
            [U_xx2(k) U_xy2(k) U_xz2(k) U_yy2(k) U_yz2(k) U_zz2(k)] = GetDerivative_6dim(mu, x2(:,k));
            [U_xx3(k) U_xy3(k) U_xz3(k) U_yy3(k) U_yz3(k) U_zz3(k)] = GetDerivative_6dim(mu, x3(:,k));
            [U_xx4(k) U_xy4(k) U_xz4(k) U_yy4(k) U_yz4(k) U_zz4(k)] = GetDerivative_6dim(mu, x4(:,k));
            [U_xx5(k) U_xy5(k) U_xz5(k) U_yy5(k) U_yz5(k) U_zz5(k)] = GetDerivative_6dim(mu, x5(:,k));
            [U_xx6(k) U_xy6(k) U_xz6(k) U_yy6(k) U_yz6(k) U_zz6(k)] = GetDerivative_6dim(mu, x6(:,k));
            [U_xx7(k) U_xy7(k) U_xz7(k) U_yy7(k) U_yz7(k) U_zz7(k)] = GetDerivative_6dim(mu, x7(:,k));
            
            dB1 = zeros(10*7, 8);
            dB1(1,1) = 1;
            dB1(4:6, 5) = [Delta_t(k)*U_xx1(k); Delta_t(k)*U_xy1(k); Delta_t(k)*U_xz1(k)]/2;
            dB1(7+2, 1) = 1;
            dB1(7+4:7+6, 5) = [Delta_t(k)*U_xy1(k); Delta_t(k)*U_yy1(k); Delta_t(k)*U_yz1(k)]/2;
            dB1(14+3, 1) = 1;
            dB1(14+4:14+6, 5) = [Delta_t(k)*U_xz1(k); Delta_t(k)*U_yz1(k); Delta_t(k)*U_zz1(k)]/2;
            dB1(21+4, 1) = 1;
            dB1(21+1, 5) = Delta_t(k)*1/2;
            dB1(21+5, 5) = -Delta_t(k)*2/2;
            dB1(28+5, 1) = 1;
            dB1(28+2, 5) = Delta_t(k)*1/2;
            dB1(28+4, 5) = Delta_t(k)*2/2;
            dB1(35+6, 1) = 1;
            dB1(35+3, 5) = Delta_t(k)*1/2;
            dB1(42+7, 1) = 1;
            dB1(42+4:42+6, 5) = -Delta_t(k)/2*u1(1,k)/x1(7,k)^2*C1*uhat_inert;
            
            
            dBu1 = zeros(7*3, 8);
            dBu1(4:6, 5) = Delta_t(k)/2/x1(7,k)*C1*uhat_inert;
            dBu1(7, 5) = Delta_t(k)/2*(-1/Isp_nd/g0_nd);
            dBu1(7+4:7+6, 5) = Delta_t(k)/2*u1(1,k)/x1(7,k)*C1*[-sin(u1(2,k))*cos(u1(3,k)); cos(u1(2,k))*cos(u1(3,k)); 0];
            dBu1(14+4:14+6, 5) = Delta_t(k)/2*u1(1,k)/x1(7,k)*C1*[-cos(u1(2,k))*sin(u1(3,k)); -sin(u1(2,k))*sin(u1(3,k)); cos(u1(3,k))];
                  
            dB1(50:70, 1:8) = dBu1;
            
            dB3 = zeros(10*7, 8);
            dB3(1,2) = 1;
            dB3(4:6, 6) = [Delta_t(k)*U_xx3(k); Delta_t(k)*U_xy3(k); Delta_t(k)*U_xz3(k)]/2;
            dB3(7+2, 2) = 1;
            dB3(7+4:7+6, 6) = [Delta_t(k)*U_xy3(k); Delta_t(k)*U_yy3(k); Delta_t(k)*U_yz3(k)]/2;
            dB3(14+3, 2) = 1;
            dB3(14+4:14+6, 6) = [Delta_t(k)*U_xz3(k); Delta_t(k)*U_yz3(k); Delta_t(k)*U_zz3(k)]/2;
            dB3(21+4, 2) = 1;
            dB3(21+1, 6) = Delta_t(k)*1/2;
            dB3(21+5, 6) = -Delta_t(k)*2/2;
            dB3(28+5, 2) = 1;
            dB3(28+2, 6) = Delta_t(k)*1/2;
            dB3(28+4, 6) = Delta_t(k)*2/2;
            dB3(35+6, 2) = 1;
            dB3(35+3, 6) = Delta_t(k)*1/2;
            dB3(42+7, 2) = 1;
            dB3(42+4:42+6, 6) = -Delta_t(k)/2*u3(1,k)/x3(7,k)^2*C3*uhat_inert;
            
            dBu3 = zeros(7*3, 8);
            dBu3(4:6, 6) = Delta_t(k)/2/x3(7,k)*C3*uhat_inert;
            dBu3(7, 6) = Delta_t(k)/2*(-1/Isp_nd/g0_nd);
            dBu3(7+4:7+6, 6) = Delta_t(k)/2*u3(1,k)/x3(7,k)*C3*[-sin(u3(2,k))*cos(u3(3,k)); cos(u3(2,k))*cos(u3(3,k)); 0];
            dBu3(14+4:14+6,6) = Delta_t(k)/2*u3(1,k)/x3(7,k)*C3*[-cos(u3(2,k))*sin(u3(3,k)); -sin(u3(2,k))*sin(u3(3,k)); cos(u3(3,k))];
   
            dB3(50:70, 1:8) = dBu3;
            
            dB5 = zeros(10*7, 8);
            dB5(1,3) = 1;
            dB5(4:6, 7) = [Delta_t(k)*U_xx5(k); Delta_t(k)*U_xy5(k); Delta_t(k)*U_xz5(k)]/2;
            dB5(7+2, 3) = 1;
            dB5(7+4:7+6, 7) = [Delta_t(k)*U_xy5(k); Delta_t(k)*U_yy5(k); Delta_t(k)*U_yz5(k)]/2;
            dB5(14+3, 3) = 1;
            dB5(14+4:14+6, 7) = [Delta_t(k)*U_xz5(k); Delta_t(k)*U_yz5(k); Delta_t(k)*U_zz5(k)]/2;
            dB5(21+4, 3) = 1;
            dB5(21+1, 7) = Delta_t(k)*1/2;
            dB5(21+5, 7) = -Delta_t(k)*2/2;
            dB5(28+5, 3) = 1;
            dB5(28+2, 7) = Delta_t(k)*1/2;
            dB5(28+4, 7) = Delta_t(k)*2/2;
            dB5(35+6, 3) = 1;
            dB5(35+3, 7) = Delta_t(k)*1/2;
            dB5(42+7, 3) = 1;
            dB5(42+4:42+6, 7) = -Delta_t(k)/2*u5(1,k)/x5(7,k)^2*C5*uhat_inert;
            
            
            dBu5 = zeros(7*3, 8);
            dBu5(4:6, 7) = Delta_t(k)/2/x5(7,k)*C5*uhat_inert;
            dBu5(7, 7) = Delta_t(k)/2*(-1/Isp_nd/g0_nd);
            dBu5(7+4:7+6, 7) = Delta_t(k)/2*u5(1,k)/x5(7,k)*C5*[-sin(u5(2,k))*cos(u5(3,k)); cos(u5(2,k))*cos(u5(3,k)); 0];
            dBu5(14+4:14+6, 7) = Delta_t(k)/2*u5(1,k)/x5(7,k)*C5*[-cos(u5(2,k))*sin(u5(3,k)); -sin(u5(2,k))*sin(u5(3,k)); cos(u5(3,k))];
%    
            dB5(50:70, 1:8) = dBu5;
            
            dB7 = zeros(10*7, 8);
            dB7(1,4) = 1;
            dB7(4:6, 8) = [Delta_t(k)*U_xx7(k); Delta_t(k)*U_xy7(k); Delta_t(k)*U_xz7(k)]/2;
            dB7(7+2, 4) = 1;
            dB7(7+4:7+6, 8) = [Delta_t(k)*U_xy7(k); Delta_t(k)*U_yy7(k); Delta_t(k)*U_yz7(k)]/2;
            dB7(14+3, 4) = 1;
            dB7(14+4:14+6, 8) = [Delta_t(k)*U_xz7(k); Delta_t(k)*U_yz7(k); Delta_t(k)*U_zz7(k)]/2;
            dB7(21+4, 4) = 1;
            dB7(21+1, 8) = Delta_t(k)*1/2;
            dB7(21+5, 8) = -Delta_t(k)*2/2;
            dB7(28+5, 4) = 1;
            dB7(28+2, 8) = Delta_t(k)*1/2;
            dB7(28+4, 8) = Delta_t(k)*2/2;
            dB7(35+6, 4) = 1;
            dB7(35+3, 8) = Delta_t(k)*1/2;
            dB7(42+7, 4) = 1;
            dB7(42+4:42+6, 8) = -Delta_t(k)/2*u7(1,k)/x7(7,k)^2*C7*uhat_inert;
            
            dBu7 = zeros(7*3, 8);
            dBu7(4:6, 8) = Delta_t(k)/2/x7(7,k)*C7*uhat_inert;
            dBu7(7, 8) = Delta_t(k)/2*(-1/Isp_nd/g0_nd);
            dBu7(7+4:7+6, 8) = Delta_t(k)/2*u7(1,k)/x7(7,k)*C7*[-sin(u7(2,k))*cos(u7(3,k)); cos(u7(2,k))*cos(u7(3,k)); 0];
            dBu7(14+4:14+6, 8) = Delta_t(k)/2*u7(1,k)/x7(7,k)*C7*[-cos(u7(2,k))*sin(u7(3,k)); -sin(u7(2,k))*sin(u7(3,k)); cos(u7(3,k))];
   
            dB7(50:70, 1:8) = dBu7;
                         
                      
            df2 = zeros(7, 7);
            df2(1:3, 4:6) = eye(3);
            df2(4:6, 1:3) = [U_xx2(k), U_xy2(k), U_xz2(k);
                                  U_xy2(k), U_yy2(k), U_yz2(k);
                                  U_xz2(k), U_yz2(k), U_zz2(k)];
            df2(4:6, 4:6) = [0, 2, 0;
                                -2, 0, 0;
                                0, 0, 0];
            df2(4:6, 7) = -u2(1,k)/x2(7,k)^2*C2*uhat_inert;

            
            df4 = zeros(7, 7);
            df4(1:3, 4:6) = eye(3);
            df4(4:6, 1:3) = [U_xx4(k), U_xy4(k), U_xz4(k);
                                  U_xy4(k), U_yy4(k), U_yz4(k);
                                  U_xz4(k), U_yz4(k), U_zz4(k)];
            df4(4:6, 4:6) = [0, 2, 0;
                                -2, 0, 0;
                                0, 0, 0];
            df4(4:6, 7) = -u4(1,k)/x4(7,k)^2*C4*uhat_inert;

            df6 = zeros(7, 7);
            df6(1:3, 4:6) = eye(3);
            df6(4:6, 1:3) = [U_xx6(k), U_xy6(k), U_xz6(k);
                                  U_xy6(k), U_yy6(k), U_yz6(k);
                                  U_xz6(k), U_yz6(k), U_zz6(k)];
            df6(4:6, 4:6) = [0, 2, 0;
                                -2, 0, 0;
                                0, 0, 0];
            df6(4:6, 7) = -u6(1,k)/x6(7,k)^2*C6*uhat_inert;
            
            %% dfdx
            
            dfdx2 = zeros(49,49);

            dfdx2(1:7, 1:7) = df2;
            dfdx2(8:14, 8:14) = df2;
            dfdx2(15:21, 15:21) = df2;
            dfdx2(22:28, 22:28) = df2;
            dfdx2(29:35, 29:35) = df2;
            dfdx2(36:42, 36:42) = df2;
            dfdx2(43:49, 43:49) = df2;
            
            dfdxu2 = zeros(21,21);
            
            dfdxu2(1:7, 1:7) = df2;
            dfdxu2(8:14, 8:14) = df2;
            dfdxu2(15:21, 15:21) = df2;
            
            dfdx4 = zeros(49,49);
            
            dfdx4(1:7, 1:7) = df4;
            dfdx4(8:14, 8:14) = df4;
            dfdx4(15:21, 15:21) = df4;
            dfdx4(22:28, 22:28) = df4;
            dfdx4(29:35, 29:35) = df4;
            dfdx4(36:42, 36:42) = df4;
            dfdx4(43:49, 43:49) = df4;
            
            dfdxu4 = zeros(21,21);
            dfdxu4(1:7, 1:7) = df4;
            dfdxu4(8:14, 8:14) = df4;
            dfdxu4(15:21, 15:21) = df4;
            
            dfdx6 = zeros(49,49);
            
            dfdx6(1:7, 1:7) = df6;
            dfdx6(8:14, 8:14) = df6;
            dfdx6(15:21, 15:21) = df6;
            dfdx6(22:28, 22:28) = df6;
            dfdx6(29:35, 29:35) = df6;
            dfdx6(36:42, 36:42) = df6;
            dfdx6(43:49, 43:49) = df6;
            
            dfdxu6 = zeros(21,21);
            dfdxu6(1:7, 1:7) = df6;
            dfdxu6(8:14, 8:14) = df6;
            dfdxu6(15:21, 15:21) = df6;
            
            %% dfdu

            df2u = zeros(7, 3);
            
            df2u(4:6, 1) = 1/x2(7,k)*C2*uhat_inert;
            df2u(4:6, 2) = u2(1,k)/x2(7,k)*C2*[-sin(u2(2,k))*cos(u2(3,k)); cos(u2(2,k))*cos(u2(3,k)); 0];
            df2u(4:6, 3) = u2(1,k)/x2(7,k)*C2*[-cos(u2(2,k))*sin(u2(3,k)); -sin(u2(2,k))*sin(u2(3,k)); cos(u2(3,k))];
            df2u(7, 1) = (-1/Isp_nd/g0_nd);
            
            df4u = zeros(7, 3);
            
            df4u(4:6, 1) = 1/x4(7,k)*C4*uhat_inert;
            df4u(4:6, 2) = u4(1,k)/x4(7,k)*C4*[-sin(u4(2,k))*cos(u4(3,k)); cos(u4(2,k))*cos(u4(3,k)); 0];
            df4u(4:6, 3) = u4(1,k)/x4(7,k)*C4*[-cos(u4(2,k))*sin(u4(3,k)); -sin(u4(2,k))*sin(u4(3,k)); cos(u4(3,k))];
            df4u(7, 1) = (-1/Isp_nd/g0_nd);
            
            df6u = zeros(7, 3);
            
            df6u(4:6, 1) = 1/x6(7,k)*C6*uhat_inert;
            df6u(4:6, 2) = u6(1,k)/x6(7,k)*C6*[-sin(u6(2,k))*cos(u6(3,k)); cos(u6(2,k))*cos(u6(3,k)); 0];
            df6u(4:6, 3) = u6(1,k)/x6(7,k)*C6*[-cos(u6(2,k))*sin(u6(3,k)); -sin(u6(2,k))*sin(u6(3,k)); cos(u6(3,k))];
            df6u(7, 1) = (-1/Isp_nd/g0_nd);
            
            dfdu2 = zeros(21, 9);
            for i = 1:3
                dfdu2(7*(i-1)+1:7*i, 3*(i-1)+1:3*i) = df2u;
            end
            
            dfdu4 = zeros(21, 9);
            for i = 1:3
                dfdu4(7*(i-1)+1:7*i, 3*(i-1)+1:3*i) = df4u;
            end    
  
            dfdu6 = zeros(21, 9);
            for i = 1:3
                dfdu6(7*(i-1)+1:7*i, 3*(i-1)+1:3*i) = df6u;
            end        
 
            % Del2
            dF(21*(k-1)+1:21*(k-1)+7, 28*(k-1)+1-7*cnt:28*(k-1)+7-7*cnt) = reshape(dB1(1:49, :)*Phi_prime(1,:)'  - Delta_t(k)/2*[dfdx2*dB1(1:49, :)*Phi(1,:)'], [7,7]);
            dF(21*(k-1)+1:21*(k-1)+7, 28*(k-1)+8-7*cnt:28*(k-1)+14-7*cnt) = reshape(dB3(1:49, :)*Phi_prime(1,:)'  - Delta_t(k)/2*[dfdx2*dB3(1:49, :)*Phi(1,:)'], [7,7]);
            dF(21*(k-1)+1:21*(k-1)+7, 28*(k-1)+15-7*cnt:28*(k-1)+21-7*cnt) = reshape(dB5(1:49, :)*Phi_prime(1,:)'  - Delta_t(k)/2*[dfdx2*dB5(1:49, :)*Phi(1,:)'], [7,7]);
            dF(21*(k-1)+1:21*(k-1)+7, 28*(k-1)+22-7*cnt:28*(k-1)+28-7*cnt) = reshape(dB7(1:49, :)*Phi_prime(1,:)'  - Delta_t(k)/2*[dfdx2*dB7(1:49, :)*Phi(1,:)'], [7,7]);

            % Del4
            dF(21*(k-1)+8:21*(k-1)+14, 28*(k-1)+1-7*cnt:28*(k-1)+7-7*cnt)  = reshape(dB1(1:49, :)*Phi_prime(2,:)'  - Delta_t(k)/2*[dfdx4*dB1(1:49, :)*Phi(2,:)'], [7,7]);
            dF(21*(k-1)+8:21*(k-1)+14, 28*(k-1)+8-7*cnt:28*(k-1)+14-7*cnt)   = reshape(dB3(1:49, :)*Phi_prime(2,:)'  - Delta_t(k)/2*[dfdx4*dB3(1:49, :)*Phi(2,:)'], [7,7]);
            dF(21*(k-1)+8:21*(k-1)+14, 28*(k-1)+15-7*cnt:28*(k-1)+21-7*cnt)  = reshape(dB5(1:49, :)*Phi_prime(2,:)'  - Delta_t(k)/2*[dfdx4*dB5(1:49, :)*Phi(2,:)'], [7,7]);
            dF(21*(k-1)+8:21*(k-1)+14, 28*(k-1)+22-7*cnt:28*(k-1)+28-7*cnt) = reshape(dB7(1:49, :)*Phi_prime(2,:)'  - Delta_t(k)/2*[dfdx4*dB7(1:49, :)*Phi(2,:)'], [7,7]);
            
            % Del6
            dF(21*(k-1)+15:21*(k-1)+21, 28*(k-1)+1-7*cnt:28*(k-1)+7-7*cnt)  = reshape(dB1(1:49, :)*Phi_prime(3,:)'  - Delta_t(k)/2*[dfdx6*dB1(1:49, :)*Phi(3,:)'], [7,7]);
            dF(21*(k-1)+15:21*(k-1)+21, 28*(k-1)+8-7*cnt:28*(k-1)+14-7*cnt)  = reshape(dB3(1:49, :)*Phi_prime(3,:)'  - Delta_t(k)/2*[dfdx6*dB3(1:49, :)*Phi(3,:)'], [7,7]);
            dF(21*(k-1)+15:21*(k-1)+21, 28*(k-1)+15-7*cnt:28*(k-1)+21-7*cnt)   = reshape(dB5(1:49, :)*Phi_prime(3,:)'  - Delta_t(k)/2*[dfdx6*dB5(1:49, :)*Phi(3,:)'], [7,7]);
            dF(21*(k-1)+15:21*(k-1)+21, 28*(k-1)+22-7*cnt:28*(k-1)+28-7*cnt) = reshape(dB7(1:49, :)*Phi_prime(3,:)'  - Delta_t(k)/2*[dfdx6*dB7(1:49, :)*Phi(3,:)'], [7,7]);
            
            % From Controls
            dF(21*(k-1)+1:21*(k-1)+7, 7*(n+1) + 3*(k-1)+1:7*(n+1) + 3*k) = reshape((dB1(50:end, :) + dB3(50:end, :) + dB5(50:end, :) + dB7(50:end, :))*Phi_prime(1,:)' - Delta_t(k)/2*[dfdxu2*dB1(50:end, :)*Phi(1,:)' + dfdxu2*dB3(50:end, :)*Phi(1,:)' + dfdxu2*dB5(50:end, :)*Phi(1,:)' + dfdxu2*dB7(50:end, :)*Phi(1,:)' + dfdu2*dudu], [7, 3]);
            dF(21*(k-1)+8:21*(k-1)+14, 7*(n+1) + 3*(k-1)+1:7*(n+1) + 3*k) = reshape((dB1(50:end, :) + dB3(50:end, :) + dB5(50:end, :) + dB7(50:end, :))*Phi_prime(2,:)' - Delta_t(k)/2*[dfdxu4*dB1(50:end, :)*Phi(2,:)' + dfdxu4*dB3(50:end, :)*Phi(2,:)' + dfdxu4*dB5(50:end, :)*Phi(2,:)' + dfdxu4*dB7(50:end, :)*Phi(2,:)' + dfdu4*dudu], [7, 3]);
            dF(21*(k-1)+15:21*(k-1)+21, 7*(n+1) + 3*(k-1)+1:7*(n+1) + 3*k) = reshape((dB1(50:end, :) + dB3(50:end, :) + dB5(50:end, :) + dB7(50:end, :))*Phi_prime(3,:)'- Delta_t(k)/2*[dfdxu6*dB1(50:end, :)*Phi(3,:)' + dfdxu6*dB3(50:end, :)*Phi(3,:)' + dfdxu6*dB5(50:end, :)*Phi(3,:)' + dfdxu6*dB7(50:end, :)*Phi(3,:)' + dfdu6*dudu], [7, 3]);
            
            cnt = cnt + 1;
        end
        
        
        dF(21*s + 1:21*s+6, 1:6) = -eye(6);
        dF(21*s + 7:21*s+12, 21*s+1:21*s+6) = -eye(6);
        dF = dF';
    end




end