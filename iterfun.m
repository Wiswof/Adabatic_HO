function [Tag_inner]=iterfun(State_coll,N,FV,hbar,m,dx,n_rela,N_state,Dn_nc)

State_coll_new=Eigen_state_value(N,FV(n_rela+Dn_nc),hbar,m,dx,N_state);
%Ground_try=State_coll_new(:,1);


for ib=1:N_state
    base_inner(ib)=State_coll_new(:,ib)'*State_coll(:,ib);
end


% fprintf(['Sign of state:',mat2str(base_inner),'\n'])


State_coll_new(:,find(base_inner<0))=-State_coll_new(:,find(base_inner<0));

%%%%%%%%%%
for ib=1:N_state
    base_inner(ib)=State_coll_new(:,ib)'*State_coll(:,ib);
end



Ground_t=State_coll(:,1);
Ground_try=State_coll_new(:,1);

% 
% Tag_inner=abs((Ground_t-Ground_try)'*Ground_t);

Tag_inner=abs(1-Ground_try'*Ground_t);

% 
% Returnb=Ground_try'*Ground_t;
% fprintf(['Sign of state:',num2str(Returnb),'\n'])

%%%%%%%%%%%%%
end