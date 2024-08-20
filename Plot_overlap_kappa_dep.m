clear;clc;close all
subplot(2,3,1)
load Data0315/data_kappa=1M.mat
Nial=10^5;
Grap.fontsize=20;
Grap.numfontsize=12;
YLIM=[10^-5,10^1]*10^6

XX=(ial-1)*Nial;
semilogy(XX,abs(inner_another_t(3:2:end,:))','LineWidth', 1.2)

ax = gca;
ax.XAxis.FontSize = Grap.numfontsize;
ax.YAxis.FontSize = Grap.numfontsize;
xlim([XX(1),XX(end)])
xlabel('$ (\frac{n}{n_c}-1) \times 10^5$','Interpreter','latex','FontSize',Grap.fontsize);
%title(['Dphi_ thre=',num2str(Dphi_thre)],'location','best')
ylabel(['$|\frac{d}{d\frac{n}{n_c}}\langle\psi_0(n)|\psi_i(n)\rangle|$'],'Interpreter','latex','FontSize',Grap.fontsize);
legend(legendInfo(3:2:end),'FontSize',Grap.numfontsize,'location','best');
grid on
ylim(YLIM)
title(['$ \kappa=',num2str(Vpar.kappa/10^6/2/pi),' $[MHz]'],'Interpreter','latex','FontSize',Grap.fontsize)
%%
subplot(2,3,2)
load Data0315/data_kappa=4M.mat
Nial=10^5;
Grap.fontsize=20;
Grap.numfontsize=12;


XX=(ial-1)*Nial;
semilogy(XX,abs(inner_another_t(3:2:end,:))','LineWidth', 1.2)

ax = gca;
ax.XAxis.FontSize = Grap.numfontsize;
ax.YAxis.FontSize = Grap.numfontsize;
xlim([XX(1),XX(end)])
xlabel('$ (\frac{n}{n_c}-1) \times 10^5$','Interpreter','latex','FontSize',Grap.fontsize);
%title(['Dphi_ thre=',num2str(Dphi_thre)],'location','best')
ylabel(['$|\frac{d}{d\frac{n}{n_c}}\langle\psi_0(n)|\psi_i(n)\rangle|$'],'Interpreter','latex','FontSize',Grap.fontsize);
legend(legendInfo(3:2:end),'FontSize',Grap.numfontsize,'location','best');
grid on
ylim(YLIM)
title(['$ \kappa=',num2str(Vpar.kappa/10^6/2/pi),' $[MHz]'],'Interpreter','latex','FontSize',Grap.fontsize)
%%
subplot(2,3,4)
load Data0315/data_kappa=7M.mat
Nial=10^5;
Grap.fontsize=20;
Grap.numfontsize=12;


XX=(ial-1)*Nial;
semilogy(XX,abs(inner_another_t(3:2:end,:))','LineWidth', 1.2)

ax = gca;
ax.XAxis.FontSize = Grap.numfontsize;
ax.YAxis.FontSize = Grap.numfontsize;
xlim([XX(1),XX(end)])
xlabel('$ (\frac{n}{n_c}-1) \times 10^5$','Interpreter','latex','FontSize',Grap.fontsize);
%title(['Dphi_ thre=',num2str(Dphi_thre)],'location','best')
ylabel(['$|\frac{d}{d\frac{n}{n_c}}\langle\psi_0(n)|\psi_i(n)\rangle|$'],'Interpreter','latex','FontSize',Grap.fontsize);
legend(legendInfo(3:2:end),'FontSize',Grap.numfontsize,'location','best');
grid on
ylim(YLIM)
title(['$ \kappa=',num2str(Vpar.kappa/10^6/2/pi),' $[MHz]'],'Interpreter','latex','FontSize',Grap.fontsize)
%%
subplot(2,3,5)
load Data0315/data_kappa=10M.mat
Nial=10^5;
Grap.fontsize=20;
Grap.numfontsize=12;


XX=(ial-1)*Nial;
semilogy(XX,abs(inner_another_t(3:2:end,:))','LineWidth', 1.2)

ax = gca;
ax.XAxis.FontSize = Grap.numfontsize;
ax.YAxis.FontSize = Grap.numfontsize;
xlim([XX(1),XX(end)])
xlabel('$ (\frac{n}{n_c}-1) \times 10^5$','Interpreter','latex','FontSize',Grap.fontsize);
%title(['Dphi_ thre=',num2str(Dphi_thre)],'location','best')
ylabel(['$|\frac{d}{d\frac{n}{n_c}}\langle\psi_0(n)|\psi_i(n)\rangle|$'],'Interpreter','latex','FontSize',Grap.fontsize);
legend(legendInfo(3:2:end),'FontSize',Grap.numfontsize,'location','best');
grid on
ylim(YLIM)
title(['$ \kappa=',num2str(Vpar.kappa/10^6/2/pi),' $[MHz]'],'Interpreter','latex','FontSize',Grap.fontsize)
%%

subplot(2,3,3)
kapparange=1:10;
for ikappa=kapparange
load(['Data0315/data_kappa=',num2str(ikappa),'M.mat'])
Maxoverlap(ikappa)=max(abs(inner_another_t(3,:)))

end

plot(kapparange,Maxoverlap,'x','markersize', 5,'LineWidth',1.2)
xlabel('$\kappa/ 2 \pi$[MHz]','Interpreter','latex','FontSize',Grap.fontsize)
ylabel(['max($|\frac{d}{d\frac{n}{n_c}}\langle\psi_0(n)|\psi_2(n)\rangle|$)'],'Interpreter','latex','FontSize',Grap.fontsize);
