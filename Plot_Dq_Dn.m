close all
figure
yyaxis left;
Grap.fontsize=30;
Grap.numfontsize=16;
plot(nr_/VPar.nc, Fdq_mindn(nr_), 'LineWidth', 2)
%legend('$\frac{dq_{min}}{dn}$','Interpreter','latex','Fontsize',30)
ax1 = gca;
ax1.XAxis.FontSize = Grap.numfontsize;
ax1.YLabel.String = '[HZ]^{-1/2}';
xlabel('$\frac{n}{n_c}$','Interpreter','latex','FontSize',Grap.fontsize);
hold on
yyaxis right;
V_height_n=simplify(subs(V_0_x,[n,x],[n,dV_0(3)]));
tmp=subs(V_height_n,[kappa,Omega,g,Delta],...
         [VPar.kappa,VPar.Omega,VPar.g,VPar.Delta]);
FV_height_n=matlabFunction(tmp);
plot(nr_/VPar.nc,FV_height_n(nr_), 'LineWidth', 2)
ax2 = gca; % 获取右侧坐标轴对象

% 设置右侧坐标轴的纵轴范围
ax2.YLim = [0, 50]; % 例如，设置范围为[0, 800]
legend(["$\frac{dq_{min}}{dn}$","$V_{height}$"],'Interpreter','latex','Fontsize',Grap.fontsize);

ax2.YLabel.String = '[HZ]';

%ax1.XAxis.FontSize = [Grap.numfontsize+5,Grap.numfontsize+5]';