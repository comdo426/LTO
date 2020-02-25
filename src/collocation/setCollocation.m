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

BAdd = nan(8,5);
BAdd(1,:) = [1,1,1,1,1];
for i = 2:8
	for j = 1:5
		BAdd(i,j) = tauAdd(j+1)^(i-1);
	end
end

invA = inv(A);

Collocation.invA = invA;
Collocation.B = B;
Collocation.D = D;
Collocation.BAdd = BAdd;
Collocation.tau = tau;
Collocation.tauAdd = tauAdd;
Collocation.tauRatio = (tau' - ones(1, 7)*tau(1))/2;

end