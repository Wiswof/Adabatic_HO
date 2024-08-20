
tic
for Fac=1:0.1:5
clearvars  -except  Fac   ;clc;
%parpool("threads")
%%%%%%%%
% Degeneracy problem can be solved by doing innerproduct of wave-func. for t+dt
% and t which should close to 1 for correct ordering
% (-1 means one '-' sign missing)
%%%%%%%%
% Constants in atomic units
L = 1500; % Region size
hbar = 1; % Reduced Planck's constant
m = 1;    % Electron mass
% Define the size of the region and the number of points for discretization
N = 10^4;   % Number of points for discretization  5*10^5 is enough for L = 5*10^3
N_state=10;   % Number of states from vacuum
dx = L / N; % Discretization step
x = linspace(-L/2, L/2, N);
%% % Define the potential function
% Vpar.Omega=2*pi*145*10^3; %Hz
% Vpar.kappa=2*pi*580*10^3; %Hz
% Vpar.Delta=2*pi*10^4; %Hz
% Vpar.g=2*pi*500; %Hz
%
% Vpar.Omega=2*pi*1; %Hz
% Vpar.kappa=2*pi*10^6; %Hz
% Vpar.Delta=2*pi*10^4; %Hz
% Vpar.g=2*pi*18.75; %Hz



    Vpar.Omega=2*pi*10; %Hz
    Vpar.kappa=2*pi*10^6*Fac; %Hz
    Vpar.Delta=2*pi*10^4; %Hz
    Vpar.g=2*pi*6.1/sqrt(Vpar.Omega/2/pi); %Hz



    %Vpar.n=17900;
    Vpar.nc=Vpar.Omega*(4*Vpar.Delta^2+Vpar.kappa^2)/(16*Vpar.g^2*Vpar.Delta);
    %Vpar.n=Vpar.nc;
    %
    % FV_alpha= @(alpha) 1/2*Vpar.Omega^2*x.^2+...
    %     alpha*Vpar.nc*(Vpar.Delta^2+Vpar.kappa^2/4)/(Vpar.kappa/2).*...
    %     (   atan( (Vpar.Delta+Vpar.g*sqrt(Vpar.Omega)*x)/(Vpar.kappa/2)   )   ...
    %     +   atan( (Vpar.Delta-Vpar.g*sqrt(Vpar.Omega)*x)/(Vpar.kappa/2)   )    );


    FV_alpha1= @(alpha) (1/2*Vpar.Omega*x.^2+...
        alpha*Vpar.nc*(Vpar.Delta^2+Vpar.kappa^2/4)/(Vpar.kappa/2).*...
        (   atan( (Vpar.Delta+Vpar.g*x)/(Vpar.kappa/2)   )   ...
        +   atan( (Vpar.Delta-Vpar.g*x)/(Vpar.kappa/2)   )    ))/Vpar.Omega    ;
    FV_alpha=@(alpha) FV_alpha1(alpha)-min(FV_alpha1(alpha));
    DHDn=((FV_alpha(1)-FV_alpha(0))/Vpar.nc);

    Nalpha=100;%
    %%
    Dalpha=0.0000005/Fac^1.3;
    ial=linspace(1-Dalpha,1+2*Dalpha,Nalpha);  %
    Dn_nc=mean(diff(ial));

    [eigenfunctions,eigenvalues]=Eigen_state_value(N,FV_alpha(ial(1)-Dn_nc),hbar,m,dx,N_state);
    Gound_state_ini=eigenfunctions(:,1) ; %N_state copies


    %%%%%%%%%%%%%%%%
    Dn_nc1=10^-7/Fac;
    %%%%%%%%%%%%
    for ia=1:Nalpha
        fprintf(['Fac=',num2str(Fac),'  Process,',num2str(ia/Nalpha),'\n']);
        toc
        V_ia=FV_alpha(ial(ia));

        if ia==1
            baseiT_1=eigenfunctions;
            EnergyiT_1=eigenvalues;
        end
        [baseiT,EnergyiT]=Eigen_state_value(N,V_ia,hbar,m,dx,N_state);
        %% flip sign
        for ib=1:N_state
            base_inner(ib)=baseiT(:,ib)'*baseiT_1(:,ib);
        end
        baseiT(:,find(base_inner<0))=-baseiT(:,find(base_inner<0));
        %%
        baseiT_plusDn1=Eigen_state_value(N,FV_alpha(ial(ia)+Dn_nc1),hbar,m,dx,N_state);
        for ib=1:N_state
            base_inner(ib)=baseiT_plusDn1(:,ib)'*baseiT(:,ib);
        end

        baseiT_plusDn1(:,find(base_inner<0))=-baseiT_plusDn1(:,find(base_inner<0));


        %%
        Ref_inner_matrix=eigenfunctions'*baseiT_1;
        Ref_inner_matrix_ground=Gound_state_ini'*baseiT_1;
        % DbaseDn_nc1=(baseiT_plusDn1-baseiT)/Dn_nc1;
        % inner_DbaseDn_nc_ground=  (baseiT(:,1))'  *DbaseDn_nc1;
        tmp=zeros(1,N_state);
        tmp(1)=1;   % inner product becomes delta_{j,1}
        inner_DbaseDn_nc_ground= (baseiT(:,1)'*baseiT_plusDn1- tmp)/Dn_nc1;
        inner_another=baseiT(:,1)'.*DHDn*baseiT*Vpar.nc./(EnergyiT-EnergyiT(1))';
        inner_another(1)=0;
        %% Data memory
        %Vmin_(ia)=min(V_ia);
        %DbaseDn_nc1_t{ia} = DbaseDn_nc1;
        Energy_spectrum(:,ia)=EnergyiT;
        inner_DbaseDn_nc_t_ground_t(:,ia)=inner_DbaseDn_nc_ground';
        inner_another_t(:,ia)=inner_another';
        baseiT_1=baseiT;
        EnergyiT_1=EnergyiT;

    end

    toc
    %%
    close all
    %%%
    Grap.fontsize=20;
    Grap.numfontsize=12;
    plot(ial,(Energy_spectrum-Energy_spectrum(1,:))'*Vpar.Omega,'LineWidth', 1.2)
    grid on
    %semilogy(ial,(Energy_spectrum-Energy_spectrum(1,:))','LineWidth', 1.2)
    ax = gca;
    ax.XAxis.FontSize = Grap.numfontsize;
    ax.YAxis.FontSize = Grap.numfontsize;

    xlabel('$\frac{n}{n_c}$','Interpreter','latex','FontSize',Grap.fontsize);
    ylabel('$(E_j(n)-E_0(n))/\hbar $[Hz]','Interpreter','latex','FontSize',Grap.fontsize);
    ylim([Energy_spectrum(1,1)-Energy_spectrum(2,1),1.05*max(max(Energy_spectrum-Energy_spectrum(1,:)))]*Vpar.Omega  )
    xlim([min(ial),max(ial)])




    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for nlotp= 1:N_state
        legendInfo{nlotp} = ['i=' num2str(nlotp-1)];
    end

    figure
    % subplot(2,1,1)
    % semilogy(ial,abs(inner_DbaseDn_nc_t_ground_t)'/Vpar.nc,'LineWidth', 1.2)
    % legend(legendInfo,'location','best');
    % xlabel('$\frac{n}{n_c}$','Interpreter','latex','FontSize',Grap.fontsize);
    % ylabel(['$|\frac{d}{dn}\langle\psi_0(n)|\psi_i(n)\rangle|$'],'Interpreter','latex','FontSize',Grap.fontsize);
    % %d\frac{n}{n_c}
    %
    % subplot(2,1,2)
    %semilogy(ial,abs(inner_DbaseDn_nc_t_ground_t(3:2:end,:))','LineWidth', 1.2)
    semilogy(ial,abs(inner_another_t(3:2:end,:))','LineWidth', 1.2)

    %hold on
    %plot(ial,-inner_DgoundDn_nc_ground/Vpar.nc,'x')
    %plot(ial,-(inner_DbaseDn_nc_t_ground_t(1,:)),'o')


    ax = gca;
    ax.XAxis.FontSize = Grap.numfontsize;
    ax.YAxis.FontSize = Grap.numfontsize;
    xlim([min(ial),max(ial)])
    xlabel('$\frac{n}{n_c}$','Interpreter','latex','FontSize',Grap.fontsize);
    %title(['Dphi_ thre=',num2str(Dphi_thre)],'location','best')
    ylabel(['$|\frac{d}{d\frac{n}{n_c}}\langle\psi_0(n)|\psi_i(n)\rangle|$'],'Interpreter','latex','FontSize',Grap.fontsize);
    legend(legendInfo(3:2:end),'FontSize',Grap.numfontsize,'location','best');
    grid on
    ylim([10^-6,10^0]*10^6)

    %%%%%%%



    % % close all
    % figure
    % AAA=2*inner_DbaseDn_nc_t_ground_t(1,:)'/Vpar.nc.*Dn_nc1_';
    % BBB=inner_DbaseDn_nc_t_ground_t(3,:)'.^2.*Vpar.nc^-2.*Dn_nc1_'.^2;
    % plot(AAA./BBB)
    %%  Plot
   

    %%%%%%%%%
    figure

    subplot(2,2,1)
    plot(ial,(Energy_spectrum-Energy_spectrum(1,:))','LineWidth', 1.2)%*Vpar.Omega
    grid on
    %semilogy(ial,(Energy_spectrum-Energy_spectrum(1,:))','LineWidth', 1.2)
    ax = gca;
    ax.XAxis.FontSize = Grap.numfontsize;
    ax.YAxis.FontSize = Grap.numfontsize;

    xlabel('$\frac{n}{n_c}$','Interpreter','latex','FontSize',Grap.fontsize);
    ylabel('$(E_j(n)-E_0(n))/(\hbar \Omega) $','Interpreter','latex','FontSize',Grap.fontsize);
    ylim([Energy_spectrum(1,1)-Energy_spectrum(2,1),1.05*max(max(Energy_spectrum-Energy_spectrum(1,:)))]  )%*Vpar.Omega
    xlim([min(ial),max(ial)])


    subplot(2,2,2)

    semilogy(ial,abs(inner_DbaseDn_nc_t_ground_t(3:2:end,:))','LineWidth', 1.2)

    %hold on
    %plot(ial,-inner_DgoundDn_nc_ground/Vpar.nc,'x')
    %plot(ial,-(inner_DbaseDn_nc_t_ground_t(1,:)),'o')


    ax = gca;
    ax.XAxis.FontSize = Grap.numfontsize;
    ax.YAxis.FontSize = Grap.numfontsize;
    xlim([min(ial),max(ial)])
    xlabel('$\frac{n}{n_c}$','Interpreter','latex','FontSize',Grap.fontsize);
    %title(['Dphi_ thre=',num2str(Dphi_thre)],'location','best')
    ylabel(['$|\frac{d}{d\frac{n}{n_c}}\langle\psi_0(n)|\psi_i(n)\rangle|$'],'Interpreter','latex','FontSize',Grap.fontsize);
    legend(legendInfo(3:2:end),'FontSize',Grap.numfontsize,'location','best');
    grid on
    %ylim([10^-6,10^-2]*10^6)

    subplot(2,2,3)

    plot(x,FV_alpha(ial(end)))


    subplot(2,2,4)
    semilogy(ial,abs(inner_another_t(3:2:end,:))','LineWidth', 1.2)

    %%
    %%%%%%%
    %% Save Data
    clear x Gound_state_ini eigenfunctions DbaseDn_nc1 baseiT_plusDn1 baseiT baseiT_1
    % 获取当前时间的字符串
    current_time_str = datestr(now, 'yyyymmddHHMMSS');
    numstr = num2str(Vpar.kappa/(2*pi)/10^6);
    % 指定文件夹
    %folder_name = '<Dpsi|psi>_Data';  % 替换为你想要保存的文件夹名称
    folder_name = 'Data0322';  % 替换为你想要保存的文件夹名称

    % 构建完整的文件路径
    file_path = fullfile(folder_name, ['data_kappa=' numstr 'M.mat']);

    % 保存数据到指定文件夹中
    save(file_path);

    disp(['Data saved at: ' file_path]);

    pause(1)
end
