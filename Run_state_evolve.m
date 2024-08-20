clear;clc;close all
%parpool("threads")
tic
%%%%%%%%
% Degeneracy problem can be solved by doing innerproduct of wave-func. for t+dt
% and t which should close to 1 for correct ordering
% (-1 means one '-' sign missing)
%%%%%%%%
% Constants in atomic units
L = 400; % Region size
hbar = 1; % Reduced Planck's constant
m = 1;    % Electron mass
% Define the size of the region and the number of points for discretization
N = 10^4;   % Number of points for discretization  5*10^5 is enough for L = 5*10^3
N_state=10;   % Number of states from vacuum
dx = L / N; % Discretization step
x = linspace(-L/2, L/2, N);
%% % Define the potential function
Vpar.Omega=2*pi*145*10^3; %Hz
Vpar.kappa=2*pi*580*10^3; %Hz
Vpar.Delta=2*pi*10^4; %Hz
Vpar.g=2*pi*500; %Hz


Vpar.nc=Vpar.Omega*(4*Vpar.Delta^2+Vpar.kappa^2)/(16*Vpar.g^2*Vpar.Delta);

FV_alpha= @(alpha) (1/2*Vpar.Omega*x.^2+...
    alpha*Vpar.nc*(Vpar.Delta^2+Vpar.kappa^2/4)/(Vpar.kappa/2).*...
    (   atan( (Vpar.Delta+Vpar.g*x)/(Vpar.kappa/2)   )   ...
    +   atan( (Vpar.Delta-Vpar.g*x)/(Vpar.kappa/2)   )    ))/Vpar.Omega    ;
DHDn_nc=((FV_alpha(1)-FV_alpha(0)));

%%
Factor=1;%10 is enough
Dalpha=0.0002*Factor;
Nalpha=500*Factor
;%
ial=linspace(1-Dalpha,1+2*Dalpha,Nalpha);  %
Dn_nc=mean(diff(ial));

[eigenfunctions,eigenvalues]=Eigen_state_value(N,FV_alpha(ial(1)-Dn_nc),hbar,m,dx,N_state);
%Gound_state_ini=eigenfunctions(:,1) ; %N_state copies
%E_ref=eigenvalues(1);
%%%%%%%%%%%%%%%%
Dn_nc1=10^-5;  %for making derivitive
%% Initial condition
T=1;  %s  



C0=zeros(N_state,1);
C0(1)=1;
Ct=zeros(N_state,1);
dt=T/(Nalpha-1);
dn_ncdt=(ial(end)-ial(1))/T;

for ia=1:Nalpha
    if mod(ia, Factor) == 0
         fprintf(['Process step1,',num2str(ia/Nalpha),'\n']);
    toc
    end

   
    V_ia=FV_alpha(ial(ia));

    if ia==1
        baseiT_1=eigenfunctions;
        EnergyiT_1=eigenvalues;
        Ct_1=C0;
    end
    [baseiT,EnergyiT]=Eigen_state_value(N,V_ia,hbar,m,dx,N_state);
%% wavefunc. at ial(ia)-Dn_nc1  ,flip sign
    for ib=1:N_state
        base_inner(ib)=baseiT(:,ib)'*baseiT_1(:,ib);
    end
    baseiT(:,find(base_inner<0))=-baseiT(:,find(base_inner<0));  %control sign of wavefunc.
%% wavefunc. at ial(ia) ,flip sign
    % baseiT_plusDn1=Eigen_state_value(N,FV_alpha(ial(ia)+Dn_nc1),hbar,m,dx,N_state);
    % for ib=1:N_state
    %     base_inner(ib)=baseiT_plusDn1(:,ib)'*baseiT(:,ib);
    % end
    % 
    % baseiT_plusDn1(:,find(base_inner<0))=-baseiT_plusDn1(:,find(base_inner<0));   %control sign of wavefunc.

    %%

    % DbaseDn_nc1=(baseiT_plusDn1-baseiT)/Dn_nc1;
    % inner_DbaseDn_nc_ground=  (baseiT(:,1))'  *DbaseDn_nc1;

    % K1= (baseiT'*baseiT_plusDn1)/Dn_nc1*dn_ncdt;
    % K1(1:N_state+1:end)=0;  % diagnoal comp to 0


    K1=baseiT'*(diag(DHDn_nc))*baseiT./(-EnergyiT+(EnergyiT)')*dn_ncdt;%./(EnergyiT-EnergyiT(1))
    K1(1:N_state+1:end)=0;
    


    K2=diag(EnergyiT-EnergyiT(1))*Vpar.Omega;
    Ct(:,ia)=Ct_1+(-i*K2+K1)*Ct_1*dt;
    Ct(:,ia)=Ct(:,ia)/sqrt(Ct(:,ia)' *Ct(:,ia));
    %% Data memory
    %Vmin_(ia)=min(V_ia);
    %DbaseDn_nc1_t{ia} = DbaseDn_nc1;
    Energy_spectrum(:,ia)=EnergyiT;   %in unit of Omega
    %inner_DbaseDn_nc_t_ground_t(:,ia)=inner_DbaseDn_nc_ground';
    baseiT_1=baseiT;
    EnergyiT_1=EnergyiT;
    %Ct_1=Ct(:,ia);

end

% 
% Grap.fontsize=20;
% Grap.numfontsize=12;
% close all
% plot(ial,(Energy_spectrum-Energy_spectrum(1,:))'*Vpar.Omega,'LineWidth', 1.2)
% grid on
% %semilogy(ial,(Energy_spectrum-Energy_spectrum(1,:))','LineWidth', 1.2)
% ax = gca;
% ax.XAxis.FontSize = Grap.numfontsize;
% ax.YAxis.FontSize = Grap.numfontsize;
% 
% xlabel('$\frac{n}{n_c}$','Interpreter','latex','FontSize',Grap.fontsize);
% ylabel('$(E_j(n)-E_0(n))/\hbar $[Hz]','Interpreter','latex','FontSize',Grap.fontsize);
% ylim([Energy_spectrum(1,1)-Energy_spectrum(2,1),1.05*max(max(Energy_spectrum-Energy_spectrum(1,:)))]*Vpar.Omega  )
% xlim([min(ial),max(ial)])
%%
subplot(1,2,1)
plot(Ct(1,:))
subplot(1,2,2)
plot(Ct(3,:))
%%
figure
plot(ial,abs(Ct(1,:)).^2)
hold on
plot(ial,abs(Ct(3,:)).^2)

plot(ial,abs(Ct(5,:)).^2)



