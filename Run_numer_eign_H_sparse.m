clear;clc;close all
%%%%%%%%%%%%%%%%%%%%%%
% need to be modified dpsi/dt is not correct, better to constrain dpsi .say 
%  dphi is 10^-5

%%%%%%%%%%%%%%%%%%%%%%

%parpool("threads")
tic
%%%%%%%%
% Degeneracy problem can be solved by doing innerproduct of wave-func. for t+dt
% and t which should close to 1 for correct ordering
% (-1 means one '-' sign missing)
%%%%%%%5
% Constants in atomic units
L = 5*10^3; % Region size
hbar = 1; % Reduced Planck's constant
m = 1;    % Electron mass
% Define the size of the region and the number of points for discretization
N = 5*10^5;   % Number of points for discretization
N_state=30;   % Number of states from vacuum
dx = L / N; % Discretization step
x = linspace(-L/2, L/2, N);
%% % Define the potential function
Vpar.Omega=2*pi*1000; %Hz
Vpar.kappa=2*pi*10^6; %Hz
Vpar.Delta=2*pi*10^4; %Hz
Vpar.g=2*pi*18.75;
Vpar.n=17900;
Vpar.nc=Vpar.Omega*(4*Vpar.Delta^2+Vpar.kappa^2)/(16*Vpar.g^2*Vpar.Delta);
%Vpar.n=Vpar.nc;


V_end=1/2*Vpar.Omega*x.^2+...
    Vpar.n*(Vpar.Delta^2+Vpar.kappa^2/4)/(Vpar.kappa/2).*...
    (   atan( (Vpar.Delta+Vpar.g*x)/(Vpar.kappa/2)   )   ...
    +   atan( (Vpar.Delta-Vpar.g*x)/(Vpar.kappa/2)   )    );
V_0=1/2*Vpar.Omega*x.^2;
V_diff=V_end-V_0;

%%%%%%%%%%%%%%%%%%%%%%   V_0+0.993565*V_diff
alpha_V=0.993572;%0.993572
V_0=V_0+alpha_V*V_diff;
V_diff=V_end-V_0;
%%%%%%%%%%%%%%%%%%%%%%

WaveThre=0.01;% trangate vaule on wavefuntion for |psi|< max(|psi|)*WaveThre
rho_0=zeros(N_state,N_state);
rho_0(1,1)=1;
dt=8*10^-3;
Speed=0.04;%0.1 %Speed of potential changing
[eigenfunctions,eigenvalues]=Eigen_state_value(N,V_0,hbar,m,dx,N_state);
NT=50;%1/Speed/dt;  %number of time
%%
for iT=1:NT
    iT
    if iT==1
        baseiT_1=eigenfunctions;
        rhoiT_1=rho_0;
        EnergyiT_1=eigenvalues;
    end

%  V_0+iT*dt*Speed*V_diff
    [baseiT,EnergyiT]=Eigen_state_value(N,V_0+iT*dt*Speed*V_diff,hbar,m,dx,N_state);



    for ib=1:N_state
        base_inner(ib)=baseiT(:,ib)'*baseiT_1(:,ib);
    end
    baseiT(:,find(base_inner<0))=-baseiT(:,find(base_inner<0));
%     for ib=1:N_state
%         base_inner1(ib)=baseiT(:,ib)'*baseiT_1(:,ib);
%     end
    DbaseDt=(baseiT-baseiT_1)/dt;
    K1=DbaseDt'*baseiT_1; %matrix on density matrix
    K2=diag(EnergyiT);

    drhodt=K1*rhoiT_1+rhoiT_1*K1'-i*K2*rhoiT_1+i*rhoiT_1*K2  ;
    rho{iT}=rhoiT_1+drhodt*dt;
    Energy_spectrum(:,iT)=EnergyiT;
    Changing_speed(:,iT)=diag(K1);



    baseiT_1=baseiT;
    EnergyiT_1=EnergyiT;
    rhoiT_1=rho{iT};

end


for irho=1:length(rho)
Vacuum_comp(:,irho)=diag(rho{irho});

end

subplot(2,1,1)
plot([1:NT]*dt*Speed*(1-alpha_V)+alpha_V,real(Vacuum_comp)')
subplot(2,1,2)
plot([1:NT]*dt*Speed*(1-alpha_V)+alpha_V,(Energy_spectrum-Energy_spectrum(1,:))')
toc
figure
plot([1:NT]*dt*Speed*(1-alpha_V)+alpha_V,Changing_speed(1,:))




%%%%%%%%
%  eigenfunctions(:,1)'*eigenfunctions(:,1)=1
%%%%%%%%
% Sort the eigenvalues and corresponding eigenfunctions in ascending order
%     [eigenvalues, idx] = sort(diag(eigenvalues));
%     eigenfunctions = eigenfunctions(:, idx);

% Plot the potential function and the first few eigenfunctions



%
%Please pay attention to eigenfunction3 and eigenfunction4, they are inverted compare with quantum mechanics textbook
%for example see:

