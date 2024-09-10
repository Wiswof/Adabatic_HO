clear

% tspan_set=4*[logspace(-2, 0, 10)];
% for iT=1:length(tspan_set)
%     iT
%     tspan=[0,tspan_set(iT)];
%     density_matrix_evolve
%     PICname=['No.',sprintf('%02d', iT),'  t_0=',num2str(tspan(2)),'  Gamma=',num2str(Gamma_m),'  dn_nc',num2str(Dalpha)];
%     saveas(gcf, ['PIC20240822/WignerF/',PICname,'.png'])
%     rhot_save={rhot{1},rhot{end}};
%     save(['PIC20240822/WignerF/',PICname,'.mat'],"rhot_save","p","q","baseiT","WignerF")
%     close all
% end


Gamma_mset=[logspace(-2, -4, 10),0];
for iG=1:length(Gamma_mset)
    iG
    Gamma_m=Gamma_mset(iG);
    density_matrix_evolve
    PICname=['No.',sprintf('%02d', iG),'  t_0=',num2str(tspan(2)),'  Gamma=',num2str(Gamma_m)];
    saveas(gcf, ['PIC20240822/WignerF/',PICname,'.png'])

    close all
end







%%

plot(ial-1,Energy_spectrum)
ax = gca;
ax.XAxis.FontSize = Grap.numfontsize;
ax.YAxis.FontSize = Grap.numfontsize;
ax.ZAxis.FontSize = Grap.numfontsize;
xlabel('$\frac{\delta n}{n_c}$','Interpreter','latex','FontSize',Grap.fontsize)
ylabel('$\frac{1}{\hbar \Omega}(E_n-V_{min})$','Interpreter','latex','FontSize',Grap.fontsize)

saveas(gcf, ['PIC20240822/WignerF/','energy_spec','.png'])
close all


plot(abs(Ci_Psi_peak_ini))
ax = gca;
ax.XAxis.FontSize = Grap.numfontsize;
ax.YAxis.FontSize = Grap.numfontsize;
ax.ZAxis.FontSize = Grap.numfontsize;
xlabel('$n$','Interpreter','latex','FontSize',Grap.fontsize)
ylabel('$\langle \psi_{int}|\psi_n(t=0)\rangle $','Interpreter','latex','FontSize',Grap.fontsize)





