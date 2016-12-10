clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DFT code for Helium atom
% - tho electrons, one nucleus
% - use of Poisson eq  to compute Vh
% - use of LDA, SLATER to compute Ex
% - use of LDA, PZ     to compute Ec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = 50;        %grid points            %50, needs to be even
b = 5;         %box length [a.u]      %5

g3= g^3;       %grid in volume
p = linspace(-b,+b,g);   
h = p(2)-p(1); %distance points grid

[X,Y,Z] = meshgrid(p,p,p);
X=X(:);Y=Y(:);Z=Z(:);
R=sqrt(X.^2+Y.^2+Z.^2);
Vext = -2./R;
e=ones(g,1);
L=spdiags([e -2*e e], -1:1, g, g)/h^2;
I = speye(g);
L3 = kron(kron(L,I),I) + kron(kron(I,L), I) + kron(kron(I,I),L);

Vtot = Vext;                      %initial guess

a = exp(-R.^2/2);                 
ncomp = -2*a / sum(a) / h^3;      %compensation charges
ncomppot = -2./R.*erf(R/sqrt(2));

%% Start SCF

for iter=1:10 %not checking convergence criterion
    disp(['----------------------------- iteration #' num2str(iter) ]);
    
    [PSI, E] = eigs(-0.5*L3+spdiags(Vtot, 0, g3, g3), 1, 'sa');

    PSI = PSI / h^(3/2); %Because we want to define the PSI such that Sum[PSI*PSI*h^3(volume element)] = 1
    n = 2*PSI.^2;        %electron density

    Vx = -(3/pi)^(1/3)*n.^(1/3); %LDA exchange

    %We want to solve Poisson equation: L3 Vh = -4*pi*n;   using conjugate gradient method
    % 
    Vh = cgs(L3, -4*pi*(n+ncomp), 1e-7, 400)-ncomppot; %conj grad sparse
    Vtot = Vext + Vh + Vx;

    T = PSI'*(-0.5*L3)*PSI * h^3;                %kinetic energy for one electron

    Eext = sum(n.*Vext)*h^3;                     %integrate External potential
    Eh = 0.5 * sum(n.*Vh)*h^3;                   %integrate Hartree potential = Eext/2

    % EXCHANGE (Slater, Gaussian: 'S', Espresso: 'sla')
    Ex = sum(-(3/4)*(3/pi)^(1/3)*n.^(4/3))*h^3;  %integrate Exchange potential

    % CORRELATION (Perdew Zunger, Gaussian: 'PL', Espresso: 'pz')
    rs=((3/4/pi)./n).^(1/3);
    %  A      B      C     D       gamma   beta1  beta2 
    c=[0.0311,-0.048,0.002,-0.0116,-0.1423,1.0529,0.3334]; %Perdew Zungler coefficients for Unpolarized corr 


    Ec=0;
    for i=1:length(rs)
        if rs(i)>=1
            ec=c(5)/(1+c(6)*sqrt(rs(i))+c(7)*rs(i));
            Ec=Ec+n(i)*ec*(1+7/6+c(6)*sqrt(rs(i))*4/3*c(7)*rs(i))/(1+c(6)*sqrt(rs(i))+c(7)*rs(i));
        else %rs(i)<1
            ec=c(1)*log(rs(i))+c(2)+c(3)*rs(i)*log(rs(i))*c(4)*rs(i);
            Ec=Ec+n(i)*c(1)*log(rs(i))+c(2)-c(1)/3+2/3*c(3)*rs(i)*log(rs(i))+(2*c(4)-c(3))*rs(i)/3;
        end
    end
    Ec=Ec*h^3;
    Etot = 2*T + Eext + Eh + Ex + Ec;

    disp(['Eigenvalue         ' num2str(E,5) ]);
    disp(['Kinetic energy     ' num2str(T,5) ]);
    disp(['External energy    ' num2str(Eext,5) ]);
    disp(['Hartree energy     ' num2str(Eh,5) ]);
    disp(['Exchange energy    ' num2str(Ex,5) ]);
    disp(['Correlation energy ' num2str(Ec,5) ]);
    disp(['Total energy       ' num2str(Etot,5) ]);
    disp(' ');

end

disp('Remember: experm solution = -2.903 [a.u.]')

% g = 80;       
% b = 5;         
% ----------------------------- iteration10
% Eigenvalue         -0.50445
% Kinetic energy      1.2944
% External energy    -6.3729
% Hartree energy      1.9465
% Exchange energy    -0.83937
% Correlation energy -0.16869
% Total energy       -2.8456
% time ca. 15 min

% g = 80;       
% b = 5;  
% ----------------------------- iteration #10
% Eigenvalue         -0.48872
% Kinetic energy      1.2051
% External energy    -6.108
% Hartree energy      1.9057
% Exchange energy    -0.81969
% Correlation energy -0.1643
% Total energy       -2.776
% time ca. 1 min
