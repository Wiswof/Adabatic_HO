clear;clc;close all
%parpool("threads")
tic
%%%%%%%%
% Degeneracy problem can be solved by doing innerproduct of wave-func. for t+dt
% and t which should close to 1 for correct ordering
% (-1 means one '-' sign missing)
%%%%%%%5
% Constants in atomic units
L = 8; % Region size
hbar = 1; % Reduced Planck's constant
m = 1;    % Electron mass
% Define the size of the region and the number of points for discretization
N = 10^3;   % Number of points for discretization  5*10^5 is enough for L = 5*10^3
N_state=8;   % Number of states from vacuum
dx = L / N; % Discretization step
x = linspace(-L/2, L/2, N);
%% % Define the potential function
% Vpar.Omega=1.5*10^5*pi; %Hz
% Vpar.kappa=2*pi*1.8*10^5; %Hz
% Vpar.Delta=2*pi*10^2; %Hz
% Vpar.g=2*pi*1875;
Vpar.Omega=2*pi*145*10^3; %Hz
Vpar.kappa=2*pi*580*10^3; %Hz
Vpar.Delta=2*pi*10^4; %Hz
Vpar.g=2*pi*500; %Hz


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

%WaveThre=0.01;% truncate vaule on wavefuntion for |psi|< max(|psi|)*WaveThre
Nalpha=3000;%
%%
Dalpha=0.001;
ial=linspace(1-Dalpha,1+2*Dalpha,Nalpha);  %
Dn_nc=mean(diff(ial));

[eigenfunctions,eigenvalues]=Eigen_state_value(N,FV_alpha(ial(1)-Dn_nc),hbar,m,dx,N_state);
%eigenfunctions=Wave_select(eigenfunctions,WaveThre,N_state);
Gound_state_ini=eigenfunctions(:,1) ; %N_state copies


%%%%%%%%%%%%%%%%
Dphi_thre=10^-8  ; % threshold for <Dphi|phi> is between Dphi_thre and Dphi_thre/23*

%%%%%%%%%%%%
for ia=1:Nalpha
    fprintf(['Process,',num2str(ia/Nalpha),'\n']);

    toc
    if ia==1
        baseiT_1=eigenfunctions;

        EnergyiT_1=eigenvalues;
    end
    [baseiT,EnergyiT]=Eigen_state_value(N,FV_alpha(ial(ia)),hbar,m,dx,N_state);


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
            fprintf(['iter step1:',num2str(Tag_iter_tmp),'  Dn_nc_down:',num2str(Dn_nc_down),...
                '  Dn_nc_up:',num2str(Dn_nc_up),'\n']) %%%%%%%%%
        end
    else
        while Tag_iter_tmp<0
            Dn_nc_up=Dn_nc_up*2;
            Tag_iter_tmp=iterfun(baseiT,N,FV_alpha,hbar,m,dx,ial(ia),N_state,Dn_nc_up)-Dphi_thre;
            fprintf(['iter step1:',num2str(Tag_iter_tmp),'  Dn_nc_down:',num2str(Dn_nc_down),...
                '  Dn_nc_up:',num2str(Dn_nc_up),'\n']) %%%%%%%%%
        end
    end

    Tag_iter=iterfun(baseiT,N,FV_alpha,hbar,m,dx,ial(ia),N_state,(Dn_nc_up+Dn_nc_down)/2);
    while abs(Tag_iter-Dphi_thre)>=Dphi_thre/2

        if Tag_iter-Dphi_thre<=-Dphi_thre/2
            Dn_nc_down=(Dn_nc_down+Dn_nc_up)/2;
            Tag_iter=iterfun(baseiT,N,FV_alpha,hbar,m,dx,ial(ia),N_state,(Dn_nc_up+Dn_nc_down)/2);
            fprintf(['iter step2:',num2str(abs(Tag_iter-Dphi_thre)),'\n']);  %%%%%%%%%
        end

        if Tag_iter-Dphi_thre>=Dphi_thre/2
            Dn_nc_up=(Dn_nc_down+Dn_nc_up)/2;
            Tag_iter=iterfun(baseiT,N,FV_alpha,hbar,m,dx,ial(ia),N_state,(Dn_nc_up+Dn_nc_down)/2);
            fprintf(['iter step2:',num2str(abs(Tag_iter-Dphi_thre)),'\n']);%%%%%%%%%
        end

    end
    %%

    Dn_nc1=(Dn_nc_up+Dn_nc_down)/2;
    

    baseiT_plusDn1=Eigen_state_value(N,FV_alpha(ial(ia)+Dn_nc1),hbar,m,dx,N_state);
    for ib=1:N_state
        base_inner(ib)=baseiT_plusDn1(:,ib)'*baseiT(:,ib);
    end

    baseiT_plusDn1(:,find(base_inner<0))=-baseiT_plusDn1(:,find(base_inner<0));




    
    %Dpsi_psi_inner=DbaseDn_nc1'*baseiT;

    %%
    Ref_inner_matrix=eigenfunctions'*baseiT_1;
    Ref_inner_matrix_ground=Gound_state_ini'*baseiT_1;
    % DbaseDn_nc1=(baseiT_plusDn1-baseiT)/Dn_nc1;
    % inner_DbaseDn_nc_ground=  (baseiT(:,1))'  *DbaseDn_nc1;
    tmp=zeros(1,N_state);
    tmp(1)=1;   % inner product becomes delta_{j,1}
    inner_DbaseDn_nc_ground= (baseiT(:,1)'*baseiT_plusDn1- tmp)/Dn_nc1;
    %% another way for ground state
    C_matrix=baseiT_plusDn1(:,1)'*baseiT;
    D_ground=-sum(C_matrix(3:2:end).^2)/2;
    
    %% Data memory

    %DbaseDn_nc1_t{ia} = DbaseDn_nc1;
    inner_DgoundDn_nc_ground(ia)=D_ground/Dn_nc1;
    Dn_nc1_(ia)=Dn_nc1;
    Energy_spectrum(:,ia)=EnergyiT;

    Ref_inner(:,ia)=diag(Ref_inner_matrix);  % inner product by ref. value
    Ref_inner_ground(:,ia)=Ref_inner_matrix_ground';
    inner_DbaseDn_nc_t_ground_t(:,ia)=inner_DbaseDn_nc_ground';
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
% subplot(2,1,1)
% semilogy(ial,abs(inner_DbaseDn_nc_t_ground_t)'/Vpar.nc,'LineWidth', 1.2)
% legend(legendInfo,'location','best');
% xlabel('$\frac{n}{n_c}$','Interpreter','latex','FontSize',Grap.fontsize);
% ylabel(['$|\frac{d}{dn}\langle\psi_0(n)|\psi_i(n)\rangle|$'],'Interpreter','latex','FontSize',Grap.fontsize);
% %d\frac{n}{n_c}
% 
% subplot(2,1,2)
semilogy(ial,abs(inner_DbaseDn_nc_t_ground_t(3:2:end,:))'/Vpar.nc,'LineWidth', 1.2)

%hold on
%plot(ial,-inner_DgoundDn_nc_ground/Vpar.nc,'x')
%plot(ial,-(inner_DbaseDn_nc_t_ground_t(1,:)),'o')


ax = gca;
ax.XAxis.FontSize = Grap.numfontsize;
ax.YAxis.FontSize = Grap.numfontsize;
xlim([min(ial),max(ial)])
xlabel('$\frac{n}{n_c}$','Interpreter','latex','FontSize',Grap.numfontsize);
%title(['Dphi_ thre=',num2str(Dphi_thre)],'location','best')
ylabel(['$|\frac{d}{dn}\langle\psi_0(n)|\psi_i(n)\rangle|$'],'Interpreter','latex','FontSize',Grap.fontsize);
legend(legendInfo(1:2:end),'FontSize',Grap.numfontsize,'location','best');

% subplot(3,1,3)
% plot(ial,-2*inner_DbaseDn_nc_t_ground_t(1,:)'/Vpar.nc.*Dn_nc1_'    ,'LineWidth', 1.2)
% hold on
% plot(ial, Dn_nc1_'.^2.*Vpar.nc^-2.*inner_DbaseDn_nc_t_ground_t(3,:)'.^2,'LineWidth', 1.2)
% plot(ial, Vpar.nc^-2*inner_DbaseDn_nc_t_ground_t(5,:)'.^2,'LineWidth', 1.2)
% plot(ial, Vpar.nc^-2*inner_DbaseDn_nc_t_ground_t(7,:)'.^2,'LineWidth', 1.2)
% legend({'$2*f_0$';'$|f_2|^2$';'$|f_4|^2$';'$|f_6|^2$'},'Interpreter','latex','FontSize',Grap.fontsize)

%%%%%%%



% % close all
% figure
% AAA=2*inner_DbaseDn_nc_t_ground_t(1,:)'/Vpar.nc.*Dn_nc1_';
% BBB=inner_DbaseDn_nc_t_ground_t(3,:)'.^2.*Vpar.nc^-2.*Dn_nc1_'.^2;
% plot(AAA./BBB)






%%%%%%%%
%% Save Data
clear x Gound_state_ini eigenfunctions DbaseDn_nc1 baseiT_plusDn1 baseiT baseiT_1
% 获取当前时间的字符串
current_time_str = datestr(now, 'yyyymmddHHMMSS');

% 指定文件夹
folder_name = '<Dpsi|psi>_Data';  % 替换为你想要保存的文件夹名称

% 构建完整的文件路径
file_path = fullfile(folder_name, ['data_' current_time_str '.mat']);

% 保存数据到指定文件夹中
save(file_path);

disp(['Data saved at: ' file_path]);
