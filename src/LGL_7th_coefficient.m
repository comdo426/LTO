%{
Last modified: 10/06/2019
Created by: Beom Park
Name: LGL_7th_coefficient_calc.m
Description: 
- Follows Zhang(2014) notation
%}

function [Phi, Phi_prime, Phi_mesh_add] = LGL_7th_coefficient

tau = [-1; -sqrt(495 +66*sqrt(15))/33; -sqrt(495 -66*sqrt(15))/33; 0; sqrt(495 -66*sqrt(15))/33; sqrt(495 +66*sqrt(15))/33; 1];

tau_mesh_add = [-1; -1+0.5*(1-sqrt(495 -66*sqrt(15))/33); -1+0.5*(1+sqrt(495 -66*sqrt(15))/33); 0; + 0.5*(1-sqrt(495 -66*sqrt(15))/33);  0.5*(1+sqrt(495 -66*sqrt(15))/33); 1];
tau_mesh_remove = [-1+2*(1-sqrt(495 -66*sqrt(15))/33); -1+2*(1+sqrt(495 -66*sqrt(15))/33)-1];

% (8, 8)
A = [1, tau(1), tau(1)^2, tau(1)^3, tau(1)^4, tau(1)^5, tau(1)^6, tau(1)^7;
    1, tau(3), tau(3)^2, tau(3)^3, tau(3)^4, tau(3)^5, tau(3)^6, tau(3)^7;
    1, tau(5), tau(5)^2, tau(5)^3, tau(5)^4, tau(5)^5, tau(5)^6, tau(5)^7;
    1, tau(7), tau(7)^2, tau(7)^3, tau(7)^4, tau(7)^5, tau(7)^6, tau(7)^7;
    0, 1, 2*tau(1), 3*tau(1)^2, 4*tau(1)^3, 5*tau(1)^4, 6*tau(1)^5, 7*tau(1)^6;
    0, 1, 2*tau(3), 3*tau(3)^2, 4*tau(3)^3, 5*tau(3)^4, 6*tau(3)^5, 7*tau(3)^6;
    0, 1, 2*tau(5), 3*tau(5)^2, 4*tau(5)^3, 5*tau(5)^4, 6*tau(5)^5, 7*tau(5)^6;
    0, 1, 2*tau(7), 3*tau(7)^2, 4*tau(7)^3, 5*tau(7)^4, 6*tau(7)^5, 7*tau(7)^6];

% (3, 8)
Xi = [1, tau(2), tau(2)^2, tau(2)^3, tau(2)^4, tau(2)^5, tau(2)^6, tau(2)^7;
    1, tau(4), tau(4)^2, tau(4)^3, tau(4)^4, tau(4)^5, tau(4)^6, tau(4)^7;
    1, tau(6), tau(6)^2, tau(6)^3, tau(6)^4, tau(6)^5, tau(6)^6, tau(6)^7];

% (3, 8)
Xi_dot = [0, 1, 2*tau(2), 3*tau(2)^2, 4*tau(2)^3, 5*tau(2)^4, 6*tau(2)^5, 7*tau(2)^6;
            0, 1, 2*tau(4), 3*tau(4)^2, 4*tau(4)^3, 5*tau(4)^4, 6*tau(4)^5, 7*tau(4)^6;
            0, 1, 2*tau(6), 3*tau(6)^2, 4*tau(6)^3, 5*tau(6)^4, 6*tau(6)^5, 7*tau(6)^6];
%         
Xi_mesh_add = [1, tau_mesh_add(2), tau_mesh_add(2)^2, tau_mesh_add(2)^3, tau_mesh_add(2)^4, tau_mesh_add(2)^5, tau_mesh_add(2)^6, tau_mesh_add(2)^7;
                        1, tau_mesh_add(3), tau_mesh_add(3)^2, tau_mesh_add(3)^3, tau_mesh_add(3)^4, tau_mesh_add(3)^5, tau_mesh_add(3)^6, tau_mesh_add(3)^7;
                        1, tau_mesh_add(4), tau_mesh_add(4)^2, tau_mesh_add(4)^3, tau_mesh_add(4)^4, tau_mesh_add(4)^5, tau_mesh_add(4)^6, tau_mesh_add(4)^7;
                        1, tau_mesh_add(5), tau_mesh_add(5)^2, tau_mesh_add(5)^3, tau_mesh_add(5)^4, tau_mesh_add(5)^5, tau_mesh_add(5)^6, tau_mesh_add(5)^7;
                        1, tau_mesh_add(6), tau_mesh_add(6)^2, tau_mesh_add(6)^3, tau_mesh_add(6)^4, tau_mesh_add(6)^5, tau_mesh_add(6)^6, tau_mesh_add(6)^7];


Phi = Xi*inv(A);
Phi_prime = Xi_dot*inv(A);
Phi_mesh_add = Xi_mesh_add*inv(A);

end




