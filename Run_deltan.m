clear  ;clc; close all

tic
% Constants in atomic units
L = 1500; % Region size
hbar = 1; % Reduced Planck's constant
m = 1;    % Electron mass
% Define the size of the region and the number of points for discretization
N = 5*10^3;%10^4   % Number of points for discretization  5*10^5 is enough for L = 5*10^3
N_int=400;  % Number of points for integration on lambda
N_state=3;   % Number of states from vacuum
dx = L / N; % Discretization step
x = linspace(-L/2, L/2, N);
%% % Define the potential function


Vpar.gamma_m=2*pi*3*10^-3; %Hz

Vpar.Omega=2*pi*1; %Hz
Vpar.kappa=2*pi*10^6; %Hz
Vpar.Delta=2*pi*10^4; %Hz
Vpar.g=2*pi*6.1/sqrt(Vpar.Omega/2/pi); %Hz



%%
Vpar.nc=Vpar.Omega*(4*Vpar.Delta^2+Vpar.kappa^2)/(16*Vpar.g^2*Vpar.Delta);
FV_alpha= @(alpha) (1/2*Vpar.Omega*x.^2+...
    alpha*Vpar.nc*(Vpar.Delta^2+Vpar.kappa^2/4)/(Vpar.kappa/2).*...
    (   atan( (Vpar.Delta+Vpar.g*x)/(Vpar.kappa/2)   )   ...
    +   atan( (Vpar.Delta-Vpar.g*x)/(Vpar.kappa/2)   )    ))/Vpar.Omega    ;
DHDn=((FV_alpha(1)-FV_alpha(0))/Vpar.nc);
%%

% <dpsi|psi>(lambda) , Energy [Hz]


Foverlap= @(lambda) overlap(1+lambda,FV_alpha,N,hbar,m,dx,N_state,Vpar,DHDn);  
%1/DeltaE already included in overlap
%%
%integral(Foverlap,0,10^-5)

t_coh=(Vpar.gamma_m+4*Vpar.Omega^2/Vpar.kappa)^-1;  %s
F_lambda=@(lambda) (custom_integral(Foverlap,0,lambda,N_int)-t_coh)/t_coh;
options = optimset('PlotFcns',{@optimplotx,@optimplotfval});
delta_lambda=fzero(F_lambda,10^-5,options);


epsilon1=(t_coh*Vpar.Omega)/delta_lambda *10^-7

figure
plot(x,FV_alpha(1+delta_lambda))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%  Plot




% %%%%%%%%%
% figure
% % Grap.fontsize=20;
% Grap.numfontsize=12;

% subplot(2,2,1)
% plot(ial,(Energy_spectrum-Energy_spectrum(1,:))','LineWidth', 1.2)%*Vpar.Omega
% grid on
% %semilogy(ial,(Energy_spectrum-Energy_spectrum(1,:))','LineWidth', 1.2)
% ax = gca;
% ax.XAxis.FontSize = Grap.numfontsize;
% ax.YAxis.FontSize = Grap.numfontsize;
% 
% xlabel('$\frac{n}{n_c}$','Interpreter','latex','FontSize',Grap.fontsize);
% ylabel('$(E_j(n)-E_0(n))/(\hbar \Omega) $','Interpreter','latex','FontSize',Grap.fontsize);
% ylim([Energy_spectrum(1,1)-Energy_spectrum(2,1),1.05*max(max(Energy_spectrum-Energy_spectrum(1,:)))]  )%*Vpar.Omega
% xlim([min(ial),max(ial)])
% 
% 
% subplot(2,2,3)
% 
% plot(x,FV_alpha(ial(end)))
% 
% 
% subplot(2,2,4)
% semilogy(ial,abs(inner_another_t(3:2:end,:))','LineWidth', 1.2)


