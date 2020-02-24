function Collocation = setCollocation
%SETCOLLOCATION - sets the collocation structure based on LGL 7th order
%
%  Syntax:
%     State = LTOMAIN(System, State, Spacecraft, Option)
%
%  Description:
%     Solves the LTO problem with LGL 7th order collocation scheme
%
%  Inputs:
%		N/A
%
%  Outputs:
%     Collocation: structure that contains the collocation information
%			tau
%
%   Author: Beom Park
%   Date: 16-Feb-2020; Last revision: 16-Feb-2020

tau = [-1; -sqrt(495 +66*sqrt(15))/33; -sqrt(495 -66*sqrt(15))/33; 0; sqrt(495 -66*sqrt(15))/33; sqrt(495 +66*sqrt(15))/33; 1];

tauAdd = [-1; -1+0.5*(1-sqrt(495 -66*sqrt(15))/33); -1+0.5*(1+sqrt(495 -66*sqrt(15))/33); 0; + 0.5*(1-sqrt(495 -66*sqrt(15))/33);  0.5*(1+sqrt(495 -66*sqrt(15))/33); 1];

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
% (5, 8)       
Xi_mesh_add = [1, tauAdd(2), tauAdd(2)^2, tauAdd(2)^3, tauAdd(2)^4, tauAdd(2)^5, tauAdd(2)^6, tauAdd(2)^7;
                        1, tauAdd(3), tauAdd(3)^2, tauAdd(3)^3, tauAdd(3)^4, tauAdd(3)^5, tauAdd(3)^6, tauAdd(3)^7;
                        1, tauAdd(4), tauAdd(4)^2, tauAdd(4)^3, tauAdd(4)^4, tauAdd(4)^5, tauAdd(4)^6, tauAdd(4)^7;
                        1, tauAdd(5), tauAdd(5)^2, tauAdd(5)^3, tauAdd(5)^4, tauAdd(5)^5, tauAdd(5)^6, tauAdd(5)^7;
                        1, tauAdd(6), tauAdd(6)^2, tauAdd(6)^3, tauAdd(6)^4, tauAdd(6)^5, tauAdd(6)^6, tauAdd(6)^7];

invA = inv(A);

phi = Xi*invA;
phiPrime = Xi_dot*invA;
phiMeshAdd = Xi_mesh_add*invA;

Collocation.phi = phi;
Collocation.phiPrime = phiPrime;
Collocation.phiMeshAdd = phiMeshAdd;
Collocation.tau = tau;
% Collocation.tauAdd = tauAdd;
Collocation.tauRatio = (tau - ones(1, 7)*tau(1))/2;

end