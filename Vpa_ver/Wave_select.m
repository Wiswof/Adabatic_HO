function [Wavenew] = Wave_select(Wave,Thres,N_state)

Wavenew=abs(Wave);
%Max_value=max(abs(Wave));
for iwave=1:N_state
P_new{iwave}=find(Wavenew(:,iwave)  >=  max(Wavenew(:,iwave))*Thres);
end

Wavenew=zeros(size(Wave,1),N_state);

for iwave=1:N_state
    Wavenew(P_new{iwave},iwave)=Wave(P_new{iwave},iwave);
    Wavenew(P_new{iwave},iwave)=Wavenew(P_new{iwave},iwave)/(Wavenew(P_new{iwave},iwave)'*Wavenew(P_new{iwave},iwave));
end



end