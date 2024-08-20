clear;clc;close all
Grap.numfontsize=12;
Grap.fontsize=20;
Vpar1.Omega0=2*pi*1; %Hz
Vpar1.kappa0=2*pi*10^6; %Hz
Vpar1.Delta0=2*pi*10^4; %Hz
Vpar1.g0=2*pi*6.1/sqrt(Vpar1.Omega0/2/pi); %Hz



Mag_fac=3;
%% fig1
folder_path1='Data0329fig1/';
folder_path='Data0329fig1/';
%folder_path='Data0401_1/';
files = dir(folder_path);
files= files(4:end);
Names={files.name};
%%%%%%%%%%%%%%%%%%%


%facDelta,fackappa,facg


%%%%%%%%%%%%%%%%%%%
for ina=1:length(Names)
    str= Names{ina};
    numbers = regexp(str, '\d+(\.\d+)?', 'match');
    paras(ina,:)=str2double(numbers);
    load([folder_path,files(ina).name])
    [Mat_overlap(ina),posimax]=max(abs(inner_another_t(3,:)));
    [~,pos_3dis]=min(abs(ial-1-Mag_fac*(ial(posimax)-1)));
    Mat_overlap1(ina)=abs(inner_another_t(3,pos_3dis));
end
kappafacrange=unique(paras(:,2));
subplot(2,2,3)
hold on
for ikappa=2:2
    filenum=find(paras(:,2)==kappafacrange(ikappa));
    [xx,xxorder]=sort(paras(filenum,2)./paras(filenum,3));
    plot(xx,Mat_overlap1(filenum(xxorder)),'LineWidth', 0.8)
    legendInfo1{ikappa} = ['\kappa_0=' num2str(kappafacrange(ikappa))];
end

skip=1;
for ikappa=5:5
    filenum=find(paras(:,2)==kappafacrange(ikappa));
    [xx,xxorder]=sort(paras(filenum,2)./paras(filenum,3));

    plot(xx(90:skip:end),Mat_overlap1(filenum(xxorder(90:skip:end))),'x','LineWidth', 1.3,'MarkerSize',7)
    legendInfo1{ikappa} = ['\kappa_0=' num2str(kappafacrange(ikappa))];

end
legendInfo2{1}=legendInfo1{2};
legendInfo2{2}=legendInfo1{5};

legend(legendInfo2,'FontSize',Grap.numfontsize,'location','best')
xlabel('$\frac{\kappa_0}{g_0}$','Interpreter','latex','FontSize',Grap.fontsize+10)
ylabel(['$|\frac{d}{d\frac{n}{n_c}}\langle\psi_0(n)|\psi_2(n)\rangle|_{asymptotic}$'],'Interpreter','latex','FontSize',Grap.fontsize+10);
ax = gca;
ax.XAxis.FontSize = Grap.numfontsize+4;
ax.YAxis.FontSize = Grap.numfontsize+4;
xlim([0,5])
%% lefttop
subplot(2,2,1)
for nlotp= 1:N_state
    legendInfo{nlotp} = ['i=' num2str(nlotp-1)];
end
icenter=find(paras(:,1)==1&paras(:,2)==1&paras(:,3)==1);
load([folder_path1,files(icenter).name])
semilogy(ial-1,abs(inner_another_t(3:2:end,:))','LineWidth', 1.2)


ax = gca;
ax.XAxis.FontSize = Grap.numfontsize+4;
ax.YAxis.FontSize = Grap.numfontsize+4;
xlim([min(ial),max(ial)]-1)
xlabel('$\frac{n}{n_c}-1$','Interpreter','latex','FontSize',Grap.fontsize);
%title(['Dphi_ thre=',num2str(Dphi_thre)],'location','best')
ylabel(['$|\frac{d}{d\frac{n}{n_c}}\langle\psi_0(n)|\psi_i(n)\rangle|$'],'Interpreter','latex','FontSize',Grap.fontsize);
legend(legendInfo(3:2:end),'FontSize',Grap.numfontsize,'location','best');
grid on
ylim([10^2,10^7])
title(['$\kappa_0=1$'],'Interpreter','latex','FontSize',Grap.fontsize+2)
%% righttop
 subplot(2,2,2)
for nlotp= 1:N_state
    legendInfo{nlotp} = ['i=' num2str(nlotp-1)];
end
load(['Data0401/','Omega=1,kappa=3M,delta=10k,g=6.1.mat'])
semilogy(ial-1,abs(inner_another_t(3:2:end,:))','LineWidth', 1.2)

ax = gca;
ax.XAxis.FontSize = Grap.numfontsize+4;
ax.YAxis.FontSize = Grap.numfontsize+4;
xlim([min(ial),max(ial)]-1)
xlabel('$\frac{n}{n_c}-1$','Interpreter','latex','FontSize',Grap.fontsize);
%title(['Dphi_ thre=',num2str(Dphi_thre)],'location','best')
ylabel(['$|\frac{d}{d\frac{n}{n_c}}\langle\psi_0(n)|\psi_i(n)\rangle|$'],'Interpreter','latex','FontSize',Grap.fontsize);
legend(legendInfo(3:2:end),'FontSize',Grap.numfontsize,'location','best');
grid on
ylim([10^2,10^7])
title(['$\kappa_0=3$'],'Interpreter','latex','FontSize',Grap.fontsize+2)

%% fig2
folder_path='Data0329fig2/';
% folder_path='Data0401_2/';

files = dir(folder_path);
files= files(4:end);
Names={files.name};
%%%%%%%%%%%%%%%%%%%
%facDelta,fackappa,facg
%%%%%%%%%%%%%%%%%%%
for ina=1:length(Names)
    str= Names{ina};
    numbers = regexp(str, '\d+(\.\d+)?', 'match');
    paras(ina,:)=str2double(numbers);
    load([folder_path,files(ina).name]);
    [Mat_overlap(ina),posimax]=max(abs(inner_another_t(3,:)));
    [~,pos_3dis]=min(abs(ial-1-Mag_fac*(ial(posimax)-1)));
    Mat_overlap1(ina)=abs(inner_another_t(3,pos_3dis));
end


kappafacrange=unique(paras(:,2));
subplot(2,2,4)
hold on
for ikappa=2
    filenum=find(paras(:,2)==kappafacrange(ikappa));
    [xx,xxorder]=sort(paras(filenum,1)./paras(filenum,2));

    plot(Vpar1.Delta0/Vpar1.kappa0*xx,Mat_overlap1(filenum(xxorder)),'LineWidth', 1)

end
skip=5;
for ikappa=length(kappafacrange)
    filenum=find(paras(:,2)==kappafacrange(ikappa));
    [xx,xxorder]=sort(paras(filenum,1)./paras(filenum,2));

    plot(Vpar1.Delta0/Vpar1.kappa0*xx(1:skip:end),Mat_overlap1(filenum(xxorder(1:skip:end)))...
        ,'x','LineWidth', 1.3,'MarkerSize',7)

end
xlim([0,0.05])
ylim([0,11*10^5])
xlabel('$\frac{\Delta}{\kappa}$','Interpreter','latex','FontSize',Grap.fontsize+10)
ylabel(['$|\frac{d}{d\frac{n}{n_c}}\langle\psi_0(n)|\psi_2(n)\rangle|_{asymptotic}$'],'Interpreter','latex','FontSize',Grap.fontsize+10);
legendInfo3{1}=['\kappa_0=1' ]
legendInfo3{2}=['\kappa_0=2.5' ]
legend(legendInfo3,'FontSize',Grap.numfontsize,'location','best')

ax = gca;
ax.XAxis.FontSize = Grap.numfontsize+4;
ax.YAxis.FontSize = Grap.numfontsize+4;
% figure;
% contour(facrange, facrange, Mat_overlap);
% %colorbar; % 显示颜色条
%
% xlabel('Fac_{kappa}');
% ylabel('Fac_{delta}');