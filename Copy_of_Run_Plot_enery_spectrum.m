clear;clc;close all
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
N_state=10;   % Number of states from vacuum
dx = L / N; % Discretization step
x = linspace(-L/2, L/2, N);
%% % Define the potential function
Vpar.Omega=1.5*10^5*pi; %Hz
Vpar.kappa=2*pi*1.8*10^5; %Hz
Vpar.Delta=2*pi*10^2; %Hz
Vpar.g=2*pi*1875;
%Vpar.n=17900;
Vpar.nc=Vpar.Omega*(4*Vpar.Delta^2+Vpar.kappa^2)/(16*Vpar.g^2*Vpar.Delta);
%Vpar.n=Vpar.nc;

FV_alpha= @(alpha) 1/2*Vpar.Omega*x.^2+...
    alpha*Vpar.nc*(Vpar.Delta^2+Vpar.kappa^2/4)/(Vpar.kappa/2).*...
    (   atan( (Vpar.Delta+Vpar.g*x)/(Vpar.kappa/2)   )   ...
    +   atan( (Vpar.Delta-Vpar.g*x)/(Vpar.kappa/2)   )    );
% V_end=1/2*Vpar.Omega*x.^2+...
%     Vpar.n*(Vpar.Delta^2+Vpar.kappa^2/4)/(Vpar.kappa/2).*...
%     (   atan( (Vpar.Delta+Vpar.g*x)/(Vpar.kappa/2)   )   ...
%     +   atan( (Vpar.Delta-Vpar.g*x)/(Vpar.kappa/2)   )    );

WaveThre=0.01;% truncate vaule on wavefuntion for |psi|< max(|psi|)*WaveThre
Nalpha=200;%
%%
Dalpha=0.001;
ial=linspace(1-Dalpha,1+2*Dalpha,Nalpha);  %
Dn_nc=mean(diff(ial));
%Dn_nc1=10^-7;% Dn_nc1 at most  10^-8
[eigenfunctions,eigenvalues]=Eigen_state_value(N,FV_alpha(ial(1)-Dn_nc),hbar,m,dx,N_state);
%eigenfunctions=Wave_select(eigenfunctions,WaveThre,N_state);
Gound_state_ini=eigenfunctions(:,1)*ones(1,N_state) ; %N_state copies
Dphi_thre=10^-4  ; % threshold for <Dphi|phi>
for ia=1:Nalpha
    ia/Nalpha
    toc
    if ia==1
        baseiT_1=eigenfunctions;

        EnergyiT_1=eigenvalues;
    end
    %  V_0+iT*dt*Speed*V_diff
    [baseiT,EnergyiT]=Eigen_state_value(N,FV_alpha(ial(ia)),hbar,m,dx,N_state);

    %baseiT=Wave_select(baseiT,WaveThre,N_state);
    %baseiT_plusDn1=Wave_select(baseiT_plusDn1,WaveThre,N_state);

    for ib=1:N_state
        base_inner(ib)=baseiT(:,ib)'*baseiT_1(:,ib);
    end

    baseiT(:,find(base_inner<0))=-baseiT(:,find(base_inner<0));
    %%



    Tag_iter=iterfun(baseiT,N,FV_alpha,hbar,m,dx,ial(ia),N_state,Dn_nc);
    Dn_nc_up=Dn_nc;
    Dn_nc_down=Dn_nc;
    Tag_iter_tmp=Tag_iter-Dphi_thre;
    if  Tag_iter-Dphi_thre>0
        while Tag_iter_tmp>0
            Dn_nc_down=Dn_nc_down/2;
            Tag_iter_tmp=iterfun(baseiT,N,FV_alpha,hbar,m,dx,ial(ia),N_state,Dn_nc_down)-Dphi_thre;
        end
    else
        while Tag_iter_tmp<0
            Dn_nc_up=Dn_nc_up*2;
            Tag_iter_tmp=iterfun(baseiT,N,FV_alpha,hbar,m,dx,ial(ia),N_state,Dn_nc_up)-Dphi_thre;
        end
    end

    Tag_iter=iterfun(baseiT,N,FV_alpha,hbar,m,dx,ial(ia),N_state,(Dn_nc_up+Dn_nc_down)/2)；
    while abs(Tag_iter-Dphi_thre)>=Dphi_thre/10

        if Tag_iter-Dphi_thre<=-Dphi_thre/10
            Dn_nc_down=(Dn_nc_down+Dn_nc_up)/2；
            Tag_iter=iterfun(baseiT,N,FV_alpha,hbar,m,dx,ial(ia),N_state,(Dn_nc_up+Dn_nc_down)/2)；

        end

        if Tag_iter-Dphi_thre>=Dphi_thre/10
            Dn_nc_up=(Dn_nc_down+Dn_nc_up)/2；
            Tag_iter=iterfun(baseiT,N,FV_alpha,hbar,m,dx,ial(ia),N_state,(Dn_nc_up+Dn_nc_down)/2)；

        end

    end
    % iter
    Dn_nc1=(Dn_nc_up+Dn_nc_down)/2;


    baseiT_plusDn1=Eigen_state_value(N,FV_alpha(ial(ia)+Dn_nc1),hbar,m,dx,N_state);

    

    DbaseDn_nc1=-(baseiT-baseiT_plusDn1)/Dn_nc1;

    %K1_n_nc=Dbase'*baseiT_1/Dn_nc1;



    Ref_inner_matrix=eigenfunctions'*baseiT_1;
    Ref_inner_matrix_ground=Gound_state_ini'*baseiT_1;

    %     for ib=1:N_state
    %         base_inner1(ib)=baseiT(:,ib)'*baseiT_1(:,ib);
    %     end

    Energy_spectrum(:,ia)=EnergyiT;
    %Rela_inner(:,ia)=diag(K1_n_nc);
    Ref_inner(:,ia)=diag(Ref_inner_matrix);  % inner product by ref. value
    Ref_inner_ground(:,ia)=diag(Ref_inner_matrix_ground);
    baseiT_1=baseiT;
    EnergyiT_1=EnergyiT;

end

toc

%%%
Grap.fontsize=20;
Grap.numfontsize=12;
close all
plot(ial,(Energy_spectrum-Energy_spectrum(1,:))','LineWidth', 1.2)
%semilogy(ial,(Energy_spectrum-Energy_spectrum(1,:))','LineWidth', 1.2)
ax = gca;
ax.XAxis.FontSize = Grap.numfontsize;
ax.YAxis.FontSize = Grap.numfontsize;

xlabel('$\frac{n}{n_c}$','Interpreter','latex','FontSize',Grap.fontsize);
ylabel('$(E_j(n)-E_0(n))/\hbar $[Hz]','Interpreter','latex','FontSize',Grap.fontsize);
ylim([Energy_spectrum(1,1)-Energy_spectrum(2,1),1.05*max(max(Energy_spectrum-Energy_spectrum(1,:)))]  )
xlim([min(ial),max(ial)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%plot(ial, Rela_inner')
plot(ial, Ref_inner','LineWidth', 1.2)

xlim([min(ial),max(ial)])

for nlotp= 1:N_state
    legendInfo{nlotp} = ['i=' num2str(nlotp-1)];
end

xlabel('$\frac{n}{n_c}$','Interpreter','latex','FontSize',Grap.fontsize);
ylabel(['$\langle\psi_i(\frac{n}{n_c}=',num2str(ial(1)-Dn_nc),')|\psi_i(n)\rangle$'],'Interpreter','latex','FontSize',Grap.fontsize);
legend(legendInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%plot(ial, Rela_inner')
plot(ial, abs(Ref_inner_ground)','LineWidth', 1.2)

xlim([min(ial),max(ial)])

for nlotp= 1:N_state
    legendInfo{nlotp} = ['i=' num2str(nlotp-1)];
end

xlabel('$\frac{n}{n_c}$','Interpreter','latex','FontSize',Grap.fontsize);
ylabel(['$|\langle\psi_0(\frac{n}{n_c}=',num2str(ial(1)-Dn_nc),')|\psi_i(n)\rangle|$'],'Interpreter','latex','FontSize',Grap.fontsize);
legend(legendInfo);






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

