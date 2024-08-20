
    clear;clc;close all
    % Constants in atomic units
    %V0 = 0.37;
    a0 = 10^2;
    hbar = 1; % Reduced Planck's constant
    m = 1;    % Electron mass
    e = 1;    % Electron charge

    % Define the size of the region and the number of points for discretization
    L = 4 * a0; % Region size
    N = 5*10^3;   % Number of points for discretization
    dx = L / N; % Discretization step

    % Define the potential function 
    x = linspace(-L/2, L/2, N);
    
    V=0.5*x.^2;
    % Construct the kinetic energy matrix T using finite difference approximation
    d2dx2 = -2*diag(ones(1, N)) + diag(ones(1, N-1), 1) + diag(ones(1, N-1), -1);
    d2dx2 = d2dx2 / dx^2;
    T = -hbar^2 / (2*m) * d2dx2;

    % Construct the potential energy matrix V
    V_matrix = diag(V);

    % Form the Hamiltonian matrix H = T + V
    H = T + V_matrix;

    % Solve the eigenvalue problem for the Hamiltonian matrix
    [eigenfunctions, eigenvalues] = eig(H);

    % Sort the eigenvalues and corresponding eigenfunctions in ascending order
    [eigenvalues, idx] = sort(diag(eigenvalues));
    eigenfunctions = eigenfunctions(:, idx);

    % Plot the potential function and the first few eigenfunctions
    num_plot_eigenfunctions = 5;
    figure;
    % Define custom colors for the eigenfunctions
    colors = {'y', 'r', 'g', 'k', 'm'};
    plot(x, V, 'b', 'LineWidth', 2);
    hold on;

    for i = 1:num_plot_eigenfunctions
        eigenfunction = eigenfunctions(:, i);
        eigenfunction = eigenfunction / max(abs(eigenfunction)); % Normalize the eigenfunction
        plot(x, eigenfunction + eigenvalues(i), 'Color', colors{i}, 'LineWidth', 1.5);
        %plot(x, eigenfunction + eigenvalues(i), 'r', 'LineWidth', 1.5);
    end
    
    hold off;
    xlabel('x (Atomic Units)');
    ylabel('Energy (Atomic Units)');
    title('Potential Function and Eigenfunctions');
    % Place the legend in the southeast corner
    legend('V(x)', 'Eigenfunction 1', 'Eigenfunction 2', 'Eigenfunction 3', 'Eigenfunction 4', 'Eigenfunction 5', 'Location', 'southeast');
    legend('boxoff');
    xlim([-4 4]);
      ylim([0 10]);
%     
%Please pay attention to eigenfunction3 and eigenfunction4, they are inverted compare with quantum mechanics textbook
%for example see:

