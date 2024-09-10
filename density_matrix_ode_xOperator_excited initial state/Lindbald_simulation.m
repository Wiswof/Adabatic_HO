clear  ;clc; close all
% Check the save folder for step_1 first



%%  K matrix
tic
L = 250; % Region size
hbar = 1; % Reduced Planck's constant
m = 1;    % Electron mass
% Define the size of the region and the number of points for discretization
N = 2*10^4;   % Number of points for discretization  5*10^5 is enough for L = 5*10^3
N_state=100;%0;   % Number of states from vacuum
Nalpha=1000;% 1000
dx = L / N; % Discretization step
x = linspace(-L/2, L/2, N);
Vpar.Omega=2*pi*1*80; %Hz
Vpar.kappa=2*pi*10^6*0.1; %Hz
Vpar.Delta=2*pi*10^4*0.01; %Hz
Vpar.g=2*pi*6.1/sqrt(Vpar.Omega/2/pi)*20; %Hz


%%


Dalpha=0.0004*2; %size of alpha region %0.002
Vpar.nc=Vpar.Omega*(4*Vpar.Delta^2+Vpar.kappa^2)/(16*Vpar.g^2*Vpar.Delta);
FV_alpha1= @(alpha) (1/2*Vpar.Omega*x.^2+...
    alpha*Vpar.nc*(Vpar.Delta^2+Vpar.kappa^2/4)/(Vpar.kappa/2).*...
    (   atan( (Vpar.Delta+Vpar.g*x)/(Vpar.kappa/2)   )   ...
    +   atan( (Vpar.Delta-Vpar.g*x)/(Vpar.kappa/2)   )    ))/Vpar.Omega    ;
FV_alpha=@(alpha) FV_alpha1(alpha)-min(FV_alpha1(alpha));
DHDn=((FV_alpha(1)-FV_alpha(0))/Vpar.nc);
DHDn=DHDn-min(DHDn);

ial=linspace(1,1+Dalpha,Nalpha);  %
Dn_nc=mean(diff(ial));

[eigenfunctions,eigenvalues]=Eigen_state_value(N,FV_alpha(ial(1)),hbar,m,dx,N_state);  %This is nomalized energy
Gound_state_ini=eigenfunctions(:,1) ; %N_state copies
inner_0_poi=find(mod([1:N_state]+[1:N_state]',2)==1);


for ia=1:Nalpha
    fprintf(['  Process,',num2str(ia/Nalpha),'\n']);
    toc
    V_ia=FV_alpha(ial(ia));

    if ia==1
        baseiT_1=eigenfunctions;
        EnergyiT_1=eigenvalues;
    end
    [baseiT,EnergyiT]=Eigen_state_value(N,V_ia,hbar,m,dx,N_state);
    baseiT=Eigenfun_Reorder(baseiT,EnergyiT);
    %% flip sign
    for ib=1:N_state
        base_inner(ib)=baseiT(:,ib)'*baseiT_1(:,ib);
    end
    baseiT(:,find(base_inner<0))=-baseiT(:,find(base_inner<0));
    Opera_q=(baseiT'*diag(x)*baseiT);
    %%

    tmp=zeros(1,N_state);
    tmp(1)=1;   % inner product becomes delta_{j,1}

    %inner_matrx=baseiT(:,1)'.*DHDn*baseiT*Vpar.nc./(EnergyiT-EnergyiT(1))';
    inner_matrix=baseiT'.*DHDn*baseiT*Vpar.nc./(EnergyiT'-EnergyiT);
    %inner_matrx(1)=0
    inner_matrix(1:N_state+1:end)=0;
    inner_matrix(inner_0_poi)=0;

    %% Data memory

    Energy_spectrum(:,ia)=EnergyiT;
    inner_matrix_t{ia}=inner_matrix;
    Opera_q_t{ia}=Opera_q;

    baseiT_1=baseiT;
    EnergyiT_1=EnergyiT;
end


save('density_matrix_ode_xOperator_excited initial state/step1.mat')

%% density matrix


%density_matrix_evolve








