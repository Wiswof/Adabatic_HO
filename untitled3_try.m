close all

ii=1
for ip=3:2:10
    subplot(2,2,ii);
plot(baseiT_plusDn1(:,ip))

legend([['N=',num2str(ip)];"Previous ground state"])
hold on
plot(baseiT(:,1),'k-.')
ii=ii+1
end
