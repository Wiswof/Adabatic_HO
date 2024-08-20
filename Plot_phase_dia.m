

Grap.legendsize=15;
Grap.fontsize=15;
Grap.labelsize=20;
Grap.numfontsize=12;

delta_height=3;
hold on
delta_=linspace(0.5,delta_height,1000);
nr=16*delta_.^2./(4*delta_.^2+1).^2;
xlim([0,1.5])
x_fill.P1 = [1, nr,0,0, fliplr(nr)];  % x坐标，fliplr用于反转nr的顺序
y_fill.P1 = [0, delta_,delta_height,0 , zeros(size(nr))];  % y坐标

x_fill.P2 = [ nr,     1   fliplr(nr)];  % x坐标，fliplr用于反转nr的顺序
y_fill.P2 = [ delta_,  delta_height,  delta_height*ones(size(nr))];  % y坐标

x_fill.P3 = [1 1.5 1.5 1];
y_fill.P3 = [0,0 delta_height delta_height];


fill(x_fill.P1, y_fill.P1, 'c', 'FaceAlpha', 0.3, 'EdgeColor', 'none');  % 'c'表示青色，'FaceAlpha'设置透明度
fill(x_fill.P2, y_fill.P2, 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');  % 'c'表示青色，'FaceAlpha'设置透明度
fill(x_fill.P3, y_fill.P3, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');  % 'c'表示青色，'FaceAlpha'设置透明度
ax = gca;
ax.XAxis.FontSize = Grap.numfontsize;
ax.YAxis.FontSize = Grap.numfontsize;
ax.ZAxis.FontSize = Grap.numfontsize;

%plot(nr,delta_, 'LineWidth', 2, 'Color', 'k')
%plot([1, 1], [0, delta_height], 'LineWidth', 2, 'Color', 'k'); % 加粗并设置颜色
text(1.1, delta_height/1.5, {'2 global minima ','at $q=\pm q_{min}$'}, 'FontSize', Grap.fontsize,'Interpreter','latex'); % 第一个相注释
text(nr(end) , delta_height - 2.5, '1 global minimum at $q=0$', 'FontSize', Grap.fontsize,'Interpreter','latex'); % 第二个相注释
text(nr(end)+0.3 , delta_height - 1, {'3 global minima','at $q=0,\pm q_{min}$'}, 'FontSize', Grap.fontsize,'Interpreter','latex'); % 第二个相注释
xlabel('$\frac{n}{n_c} $','Interpreter','latex','FontSize',Grap.labelsize);
ylabel('$\frac{\Delta}{\kappa}$','Interpreter','latex','FontSize',Grap.labelsize,'Rotation',0);
%%%%%%%%%%%%%%%%%
%%
% Vpar.Omega=2*pi; %Hz
% Vpar.kappa=2*pi*10^6; %Hz
% Vpar.Delta=2*pi*10^4; %Hz
% Vpar.g=2*pi*18.75;
% Vpar.n=17900;
% Vpar.nc=Vpar.Omega*(4*Vpar.Delta^2+Vpar.kappa^2)/(16*Vpar.g^2*Vpar.Delta);
% 
% 
% nabs=nr*Vpar.nc;
% figure
% xlim([0,1.5]*Vpar.nc)
% hold on
% x_fill.P1 = [1*Vpar.nc, nabs,0,0, fliplr(nabs)];  % x坐标，fliplr用于反转nr的顺序
% y_fill.P1 = [0, delta_,delta_height,0 , zeros(size(nabs))];  % y坐标
% 
% x_fill.P2 = [ nabs,     1*Vpar.nc   fliplr(nabs)];  % x坐标，fliplr用于反转nr的顺序
% y_fill.P2 = [ delta_,  delta_height,  delta_height*ones(size(nabs))];  % y坐标
% 
% x_fill.P3 = [1 1.5 1.5 1]*Vpar.nc;
% y_fill.P3 = [0,0 delta_height delta_height];
% 
% 
% fill(x_fill.P1, y_fill.P1, 'c', 'FaceAlpha', 0.3, 'EdgeColor', 'none');  % 'c'表示青色，'FaceAlpha'设置透明度
% fill(x_fill.P2, y_fill.P2, 'y', 'FaceAlpha', 0.3, 'EdgeColor', 'none');  % 'c'表示青色，'FaceAlpha'设置透明度
% fill(x_fill.P3, y_fill.P3, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');  % 'c'表示青色，'FaceAlpha'设置透明度
% ax = gca;
% ax.XAxis.FontSize = Grap.numfontsize;
% ax.YAxis.FontSize = Grap.numfontsize;
% ax.ZAxis.FontSize = Grap.numfontsize;
% 
% %plot(nr,delta_, 'LineWidth', 2, 'Color', 'k')
% %plot([1, 1], [0, delta_height], 'LineWidth', 2, 'Color', 'k'); % 加粗并设置颜色
% text(1.1*Vpar.nc, delta_height/1.5, {'2 global minima ','at $q=\pm q_{min}$'}, 'FontSize', Grap.fontsize,'Interpreter','latex'); % 第一个相注释
% text(nr(end)*Vpar.nc , delta_height - 2.5, '1 global minimum at $q=0$', 'FontSize', Grap.fontsize,'Interpreter','latex'); % 第二个相注释
% text((nr(end)+0.3)*Vpar.nc , delta_height - 1, {'3 global minima','at $q=0,\pm q_{min}$'}, 'FontSize', Grap.fontsize,'Interpreter','latex'); % 第二个相注释
% xlabel('$n $','Interpreter','latex','FontSize',Grap.labelsize);
% ylabel('$\frac{\Delta}{\kappa}$','Interpreter','latex','FontSize',Grap.labelsize,'Rotation',0);