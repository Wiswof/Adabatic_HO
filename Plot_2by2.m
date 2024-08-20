clear;clc;close all
subplot(2,2,1)
load(['Data0324/Figure_0324']);
semilogy(ial,abs(inner_another_t(3:2:end,:))','LineWidth', 1.2)


ax = gca;
ax.XAxis.FontSize = Grap.numfontsize+4;
ax.YAxis.FontSize = Grap.numfontsize+4;
xlim([min(ial),max(ial)])
xlabel('$\frac{n}{n_c}$','Interpreter','latex','FontSize',Grap.fontsize);
%title(['Dphi_ thre=',num2str(Dphi_thre)],'location','best')
ylabel(['$|\frac{d}{d\frac{n}{n_c}}\langle\psi_0(n)|\psi_i(n)\rangle|$'],'Interpreter','latex','FontSize',Grap.fontsize);
legend(legendInfo(3:2:end),'FontSize',Grap.numfontsize,'location','best');
grid on
ylim([10^-4,10^0]*10^7)




%% kappa
subplot(2,2,2)
clear
kapparange=1:0.1:5;
for ikappa=1:length(kapparange)
load(['Data0322/data_kappa=',num2str(kapparange(ikappa)),'M.mat'])
Maxoverlap(ikappa)=max(abs(inner_another_t(3,:)));
ikappa
end
plot(kapparange.^1.25,Maxoverlap,'x','markersize', 5,'LineWidth',1.2)
xlabel('$\kappa_0^{1.25}$','Interpreter','latex','FontSize',Grap.fontsize)
ylabel(['max($|\frac{d}{d\frac{n}{n_c}}\langle\psi_0(n)|\psi_2(n)\rangle|$)'],'Interpreter','latex','FontSize',Grap.fontsize);

ax = gca;
ax.XAxis.FontSize = Grap.numfontsize+4;
ax.YAxis.FontSize = Grap.numfontsize+4;
xlim([0,max(kapparange)])


%% Omega
subplot(2,2,3)
%subplot(2,2,4)
Omegarange=1:1:100;
%Omegarange=100:10:100;

for iOmega=Omegarange
load(['Data0322/data_Omega0=',num2str(iOmega),'.mat']);
Maxoverlap(iOmega)=max(abs(inner_another_t(3,:)));
iOmega
end

plot(Omegarange.^0.75,Maxoverlap,'x','markersize', 5,'LineWidth',1.2)
xlabel('$\Omega_0^{0.75}$','Interpreter','latex','FontSize',Grap.fontsize+10)
ylabel(['max($|\frac{d}{d\frac{n}{n_c}}\langle\psi_0(n)|\psi_2(n)\rangle|$)'],'Interpreter','latex','FontSize',Grap.fontsize+10);

ax = gca;
ax.XAxis.FontSize = Grap.numfontsize+4;
ax.YAxis.FontSize = Grap.numfontsize+4;
%xlim([0,100])
ylim([0,max(Maxoverlap)])

%% gfac
subplot(2,2,4)
clear
gfac=1:0.1:5;
%Omegarange=100:10:100;

for igfac=1:length(gfac)
load(['Data0322/data_gfac=',num2str(gfac(igfac)),'.mat']);
Maxoverlap(igfac)=max(abs(inner_another_t(3,:)));
igfac
end

plot(gfac.^-1.5,Maxoverlap,'x','markersize', 5,'LineWidth',1.2)
xlabel('$g_{fac}^{-1.5}$','Interpreter','latex','FontSize',Grap.fontsize+10)
ylabel(['max($|\frac{d}{d\frac{n}{n_c}}\langle\psi_0(n)|\psi_2(n)\rangle|$)'],'Interpreter','latex','FontSize',Grap.fontsize+10);

ax = gca;
ax.XAxis.FontSize = Grap.numfontsize+4;
ax.YAxis.FontSize = Grap.numfontsize+4;

ylim([0,max(Maxoverlap)])