


Gamma_mset=[logspace(-5, -7, 10),0];
for iG=1:length(Gamma_mset)
  Gamma_m=Gamma_mset(iG);
  density_matrix_evolve
  PICname=['No.',sprintf('%02d', iG),'  t_0=',num2str(tspan(2)),'  Gamma=',num2str(Gamma_m)];
  saveas(gcf, ['PIC20240819/WignerF/',PICname,'.png'])

  close all
end




