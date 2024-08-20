function [eigenfunctions]=Eigenfun_Reorder(eigenfunctions,eigenvalues)

for ie = find(diff(eigenvalues)==0)'

phi_sym1=eigenfunctions(:,ie)/2+flip(eigenfunctions(:,ie))/2;   
phi_asy1=eigenfunctions(:,ie)/2-flip(eigenfunctions(:,ie))/2;   
phi_sym2=eigenfunctions(:,ie+1)/2+flip(eigenfunctions(:,ie+1))/2;
phi_asy2=eigenfunctions(:,ie+1)/2-flip(eigenfunctions(:,ie+1))/2;

Norm_sym1=sqrt(phi_sym1'*phi_sym1);
Norm_sym2=sqrt(phi_sym2'*phi_sym2);


if Norm_sym1>Norm_sym2
    Norm_asy_2=sqrt(phi_asy2'*phi_asy2);
    eigenfunctions(:,ie)=phi_sym1/Norm_sym1;
    eigenfunctions(:,ie+1)=phi_asy2/Norm_asy_2;
    
else
    Norm_asy_1=sqrt(phi_asy1'*phi_asy1);
    eigenfunctions(:,ie)=phi_sym2/Norm_sym2;
    eigenfunctions(:,ie+1)=phi_asy1/Norm_asy_1;
end
end


end
