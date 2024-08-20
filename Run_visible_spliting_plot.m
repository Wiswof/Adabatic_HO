clear;clc;close all
syms x Omega kappa Delta g n epsilon real
assume(Omega>0)
assume(kappa>0)
assume(Delta>0)
assume(g>0)
assume(n>0)
V=1/2*Omega*x.^2+...
    n*(Delta^2+kappa^2/4)/(kappa/2).*...
    (   atan( (Delta+g*x)/(kappa/2)   )   ...
    +   atan( (Delta-g*x)/(kappa/2)   )    )  %(1-10^-11)

dV=diff(V,x);
%ddV=diff(dV,x);
dV_0=solve(dV,x)
n0=solve(dV_0(3)==0,n)




L = 10^4; % Region size
hbar = 1; % Reduced Planck's constant
m = 1;    % Electron mass
% Define the size of the region and the number of points for discretization
N = 5*10^5;   % Number of points for discretization
N_state=4;   % Number of states from vacuum
dx = L / N; % Discretization step
x = linspace(-L/2, L/2, N);


labelRange=20*[-1,1];


%% % Define the potential function
% Vpar.Omega=2*pi; %Hz
% Vpar.kappa=2*pi*18.75*10^2; %Hz
% Vpar.Delta=2*pi*10^-2; %Hz   %2*pi*10^6/100
% Vpar.g=2*pi*18.75;

Vpar.Omega=1.5*10^5*pi; %Hz
Vpar.kappa=2*pi*1.8*10^5; %Hz
Vpar.Delta=2*pi*10^2; %Hz
Vpar.g=2*pi*1875;
%alpha_V=0.9935695;%%17900*0.9935715 splitiing happens for alpha_V=1
dalpha=0.0001;
alpha_0=1;
Vpar.n0=double(subs(n0,[kappa,Omega,g,Delta],[Vpar.kappa,Vpar.Omega,Vpar.g,Vpar.Delta]));
for ipar=1:8
    alpha_V=alpha_0+dalpha*(ipar-1);
    Vpar.n=Vpar.n0*alpha_V;
    V=1/2*Vpar.Omega*x.^2+...
        Vpar.n*(Vpar.Delta^2+Vpar.kappa^2/4)/(Vpar.kappa/2).*...
        (   atan( (Vpar.Delta+Vpar.g*x)/(Vpar.kappa/2)   )   ...
        +   atan( (Vpar.Delta-Vpar.g*x)/(Vpar.kappa/2)   )    );  %(1-10^-11)


    d2dx2 = sparse([[1:N],[2:N],[1:N-1]],   [[1:N],[1:N-1],[2:N]],...
        [[-2*ones(1,N),1*ones(1,N-1),1*ones(1,N-1)]],N,N); %-2*diag(ones(1, N)) + diag(ones(1, N-1), 1) + diag(ones(1, N-1), -1);
    d2dx2 = d2dx2 / dx^2;
    T = -hbar^2 / (2*m) * d2dx2;

    % Construct the potential energy matrix V
    V_matrix = sparse([1:N],[1:N],V,N,N);   %diag(V);

    % Form the Hamiltonian matrix H = T + V
    H = T + V_matrix;


    tic
    % Solve the eigenvalue problem for the Hamiltonian matrix
    [eigenfunctions, eigenvalues] =  eigs(H,N_state,min(V));
    toc
    eigenvalues=diag(eigenvalues);
    %%


    
    eigenfunctions=Eigenfun_Reorder(eigenfunctions,eigenvalues);
    subplot(4,4,2*ipar-1)
    plot(x,V-min(V))
%    hold on
%    plot(0*[1:8],eigenvalues-min(V),'o')
    title(['Potential  ','i_{\alpha}=',num2str(ipar)])
    xlabel('q')
    xlim(labelRange)
    subplot(4,4,2*ipar)
    plot(x,abs(eigenfunctions(:,1)))
    title(['Ground state wavefunction  ','i_{\alpha}=',num2str(ipar)])
    xlabel('q')
    xlim(labelRange)

end
    annotation('textbox', [0 0.9 1 0.1], 'String',...
       ['$\frac{\Omega}{2\pi}$=',num2str(Vpar.Omega/2/pi),'Hz,$\frac{\Delta}{g}$=',num2str(Vpar.Delta/Vpar.g),...
        ',$\frac{\kappa}{g}$=',num2str(Vpar.kappa/Vpar.g),',$\frac{g}{2\pi}$=',num2str(Vpar.g/2/pi),...
        '$Hz,n_0$=',num2str(Vpar.n0),',$n_\alpha=(1+',num2str(dalpha),'i_\alpha)*n_0$'], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold','Interpreter','latex');
    
    
    
    
   

