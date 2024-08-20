function [Return_Value,Return_Value_overlap,Return_Value_Energy]...
            =overlap(lambda,FV_alpha,N,hbar,m,dx,N_state,Vpar,DHDn)
%lambda=n/nc
%%%%%%%%%%%%
V_ia=FV_alpha(lambda);
[baseiT,EnergyiT]=Eigen_state_value(N,V_ia,hbar,m,dx,N_state);
inner_another=baseiT(:,1)'.*DHDn*baseiT*Vpar.nc./(EnergyiT-EnergyiT(1))';
inner_another(1)=0;
Return_Value_overlap=abs(inner_another(3));
Return_Value_Energy=(EnergyiT-EnergyiT(1))*Vpar.Omega;  %Hz
Return_Value=Return_Value_overlap/Return_Value_Energy(3);
end