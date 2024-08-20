function [eigenfunctions,eigenvalues]=Eigen_state_value(N,V,hbar,m,dx,N_state)
%%
%V=0.5*x.^2;%0.5*x.^4-L^2*x.^2/10;

% Construct the kinetic energy matrix T using finite difference approximation
d2dx2 = sparse([[1:N],[2:N],[1:N-1]],   [[1:N],[1:N-1],[2:N]],...
    [[-2*ones(1,N),1*ones(1,N-1),1*ones(1,N-1)]],N,N); %-2*diag(ones(1, N)) + diag(ones(1, N-1), 1) + diag(ones(1, N-1), -1);
d2dx2 = d2dx2 / dx^2;
T = -hbar^2 / (2*m) * d2dx2;

% Construct the potential energy matrix V
V_matrix = sparse([1:N],[1:N],V,N,N);   %diag(V);

% Form the Hamiltonian matrix H = T + V
H = T + V_matrix;

% Solve the eigenvalue problem for the Hamiltonian matrix
[eigenfunctions, eigenvalues] =  eigs(H,N_state,min(V));%,'Tolerance',1e-16
eigenvalues=diag(eigenvalues);
%%
eigenfunctions=Eigenfun_Reorder(eigenfunctions,eigenvalues);




end
