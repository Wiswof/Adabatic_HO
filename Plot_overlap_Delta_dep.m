clear;clc;close all

%%

%subplot(2,2,4)
Deltafacrange=10.^[-2:0.1:1];
%Omegarange=100:10:100;

for iDeltafac=1:length(Deltafacrange)
load(['Data0325/data_facDelta=',num2str(Deltafacrange(iDeltafac)),'.mat']);
Maxoverlap(iDeltafac)=max(abs(inner_another_t(3,:)));
iDeltafac
end

semilogx(Deltafacrange*10^4,Maxoverlap,'x','markersize', 5,'LineWidth',1.2)
xlabel('$\frac{\Delta}{2\pi}$[Hz]','Interpreter','latex','FontSize',Grap.fontsize+10)
ylabel(['max($|\frac{d}{d\frac{n}{n_c}}\langle\psi_0(n)|\psi_2(n)\rangle|$)'],'Interpreter','latex','FontSize',Grap.fontsize+10);

ax = gca;
ax.XAxis.FontSize = Grap.numfontsize+6;
ax.YAxis.FontSize = Grap.numfontsize+4;
%xlim([0,100])
ylim([0,max(Maxoverlap)])