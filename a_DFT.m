clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DFT code for Hydrogen atom
% - only one electron: no Vee, no need for exchange-correlation
% - using finite difference: no basis set    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GRID
g = 50; %grid points            %50, needs to be even
b = 5;   %box length [a.u]      %5

g3= g^3; %grid in volume
p = linspace(-b,+b,g);   
h = p(2)-p(1); %distance points grid

[X, Y, Z]=meshgrid (p, p, p);
X=X(:);
Y=Y(:);
Z=Z(:);
R =sqrt(X.^2 + Y.^2 + Z.^2); %distance from origin

%% Computation
Vext = -1./R; %external potential
e = ones(g,1); %vector 
L = spdiags([e -2*e e], -1:1, g,g)/h^2; %Laplacian, extract diagonal matrice
I = speye(g); %sparse eye matrix (sparse to be faster)
L3 = kron(kron(L,I),I) + kron(kron(I,L), I) + kron(kron(I,I),L);
[PSI,E] = eigs(-0.5*L3+spdiags(Vext, 0, g3, g3),1,'sa')

PSI = PSI / h^(3/2); %Because we want to define the PSI such that Sum[PSI*PSI*h^3(volume element)] = 1

disp('Remember: exact solution = -0.5 [a.u.]')
