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

A = nan(8,8);
A(1,:) = [1,1,1,1,0,0,0,0];
for i = 2:8
	for j = 1:4
		A(i,j) = tau(2*j-1)^(i-1);
		A(i,j+4) = (i-1)*tau(2*j-1)^(i-2);
	end
end

B = nan(8,3);
B(1,:) = [1,1,1];
for i = 2:8
	for j = 1:3
		B(i,j) = tau(2*j)^(i-1);
	end
end

D = nan(8,3);
D(1,:) = [0,0,0];
for i = 2:8
	for j = 1:3
		D(i,j) = (i-1)*tau(2*j)^(i-2);		
	end
end

% (5, 8)       
Xi_mesh_add = [1, tauAdd(2), tauAdd(2)^2, tauAdd(2)^3, tauAdd(2)^4, tauAdd(2)^5, tauAdd(2)^6, tauAdd(2)^7;
                        1, tauAdd(3), tauAdd(3)^2, tauAdd(3)^3, tauAdd(3)^4, tauAdd(3)^5, tauAdd(3)^6, tauAdd(3)^7;
                        1, tauAdd(4), tauAdd(4)^2, tauAdd(4)^3, tauAdd(4)^4, tauAdd(4)^5, tauAdd(4)^6, tauAdd(4)^7;
                        1, tauAdd(5), tauAdd(5)^2, tauAdd(5)^3, tauAdd(5)^4, tauAdd(5)^5, tauAdd(5)^6, tauAdd(5)^7;
                        1, tauAdd(6), tauAdd(6)^2, tauAdd(6)^3, tauAdd(6)^4, tauAdd(6)^5, tauAdd(6)^6, tauAdd(6)^7];

invA = inv(A);

% 
% phi = Xi*invA;
% phiPrime = Xi_dot*invA;
% phiMeshAdd = Xi_mesh_add*invA;

Collocation.invA = invA;
Collocation.B = B;
Collocation.D = D;
% Collocation.phiMeshAdd = phiMeshAdd;
% Collocation.tau = tau;
% % Collocation.tauAdd = tauAdd;
Collocation.tauRatio = (tau' - ones(1, 7)*tau(1))/2;

end