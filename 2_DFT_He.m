clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DFT code for Helium atom
% - tho electrons, one nucleus
% - use of Poisson eq  (method 1) 
%       or integration (method 2) to compute Vh
% - use of LDA, SLATER to compute Ex
% - use of LDA, PZ     to compute Ec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = 20;        %grid points (must be even)
b = 4;         %box length [a.u]      

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

    %KINETIC
    T = PSI'*(-0.5*L3)*PSI * h^3;                %kinetic energy for one electron
    
    %EXTERNAL
    Eext = sum(n.*Vext)*h^3;                     %integrate External potential    

    %HARTREE
    usemethod=2;
   %(method 1) Solving Poisson equation L3*Vh=-4*pi*n using conjugate gradient method: faster!
   if usemethod==1
    Vh = cgs(L3, -4*pi*(n+ncomp), 1e-7, 400)-ncomppot; %conj grad sparse
    
   %(method 2) Numerical integration: expensive N^2 operation!
   elseif usemethod==2
    Vh=zeros(length(X),1);
    %other method
    for i=1:length(X)
        for j=1:length(X)
            if j~=i
                dX=X(i)-X(j);
                dY=Y(i)-Y(j);
                dZ=Z(i)-Z(j);
                r=sqrt(dX^2+dY^2+dZ^2);
                Vh(i)=Vh(i)+n(j)/r*h^3;
            end
        end
    end
   end
    Eh = 0.5 * sum(n.*Vh)*h^3;                   %integrate Hartree potential 
    
    %EXCHANGE (Slater, Gaussian: 'S', Espresso: 'sla')
    Vx = -(3/pi)^(1/3)*n.^(1/3); %LDA exchange
    Ex = sum(-(3/4)*(3/pi)^(1/3)*n.^(4/3))*h^3;  %integrate Exchange potential

    %CORRELATION (Perdew Zunger, Gaussian: 'PL', Espresso: 'pz')
    rs=((3/4/pi)./n).^(1/3);
    %  A      B      C     D       gamma   beta1  beta2 
    c=[0.0311,-0.048,0.002,-0.0116,-0.1423,1.0529,0.3334]; %Perdew Zungler coefficients for Unpolarized corr 

    Vc=0;
    Ec=0;
    for i=1:length(rs)
        if rs(i)>=1
            ec=c(5)/(1+c(6)*sqrt(rs(i))+c(7)*rs(i));
            Ec=Ec+ec*n(i);
            Vc=Vc+ec*(1+7/6+c(6)*sqrt(rs(i))*4/3*c(7)*rs(i))/(1+c(6)*sqrt(rs(i))+c(7)*rs(i));
        else %rs(i)<1
            ec=c(1)*log(rs(i))+c(2)+c(3)*rs(i)*log(rs(i))+c(4)*rs(i);
            Ec=Ec+ec*n(i);
            Vc=Vc+c(1)*log(rs(i))+c(2)-c(1)/3+2/3*c(3)*rs(i)*log(rs(i))+(2*c(4)-c(3))*rs(i)/3;
        end
    end
    Ec=Ec*h^3;
    
    %TOTAL
    Vtot = Vext + Vh + Vx + Vc;
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

% g = 50;        %grid points           
% b = 3;         %box length [a.u]      
% ----------------------------- iteration #10
% cgs converged at iteration 149 to a solution with relative residual 3.9e-08.
% Eigenvalue         -5182.8
% Kinetic energy     1.3425
% External energy    -6.4913
% Hartree energy     1.9964
% Exchange energy    -0.86102
% Correlation energy -0.11136
% Total energy       -2.7822
% 1min

% g = 40;        %grid points (must be even)
% b = 2;         %box length [a.u] 
% ----------------------------- iteration #10
% cgs converged at iteration 130 to a solution with relative residual 7.2e-08.
% Eigenvalue         -2871.3
% Kinetic energy     1.625
% External energy    -7.0788
% Hartree energy     2.1531
% Exchange energy    -0.96131
% Correlation energy -0.11724
% Total energy       -2.7544
% 30sec

% g = 60;        %grid points (must be even)
% b = 3;         %box length [a.u] 
% ----------------------------- iteration #10
% cgs converged at iteration 205 to a solution with relative residual 8.9e-08.
% Eigenvalue         -8942.8
% Kinetic energy     1.3625
% External energy    -6.5486
% Hartree energy     2.0053
% Exchange energy    -0.86533
% Correlation energy -0.11157
% Total energy       -2.7952
% 2min

% g = 80;        %grid points (must be even)
% b = 4;         %box length [a.u] 
% ----------------------------- iteration #10
% cgs converged at iteration 264 to a solution with relative residual 6.2e-08.
% Eigenvalue         -19767
% Kinetic energy     1.3245
% External energy    -6.4584
% Hartree energy     1.9655
% Exchange energy    -0.84767
% Correlation energy -0.11028
% Total energy       -2.8019
% 8min
