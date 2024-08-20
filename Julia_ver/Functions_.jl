

function Eigen_state_value(N, V, hbar, m, dx, N_state)
    # Construct the kinetic energy matrix T using finite difference approximation
    I_ = vcat(collect(1:N), collect(2:N), collect(1:N-1))
    J_ = vcat(collect(1:N), collect(1:N-1), collect(2:N))
    K_com = vcat(-2 * ones(BigInt, N), ones(BigInt, N - 1), ones(BigInt, N - 1))
    d2dx2 = sparse(I_, J_, K_com, N, N)
    d2dx2 = d2dx2 / dx^2
    T = -hbar^2 / (2m) * d2dx2

    # Construct the potential energy matrix V
    V_matrix = sparse(1:N, 1:N, V, N, N)

    # Form the Hamiltonian matrix H = T + V
    H = T + V_matrix
    
    println(V_matrix)

    # Solve the eigenvalue problem for the Hamiltonian matrix
    eigenvalues, eigenfunctions = eigs(H, nev=N_state, which=:SR)
    #println(eigenvalues)
    # Extract diagonal elements from eigenvalues
    eigenvalues = diagm(eigenvalues)

    # Reorder eigenfunctions
    #  eigenfunctions = Eigenfun_Reorder(eigenfunctions, eigenvalues)

    return eigenfunctions, eigenvalues
end




function Eigenfun_Reorder(eigenfunctions, eigenvalues)
    for ie in findall(diff(eigenvalues) .== 0)
        phi_sym1 = eigenfunctions[:, ie] / 2 + reverse(eigenfunctions[:, ie]) / 2
        phi_asy1 = eigenfunctions[:, ie] / 2 - reverse(eigenfunctions[:, ie]) / 2
        phi_sym2 = eigenfunctions[:, ie+1] / 2 + reverse(eigenfunctions[:, ie+1]) / 2
        phi_asy2 = eigenfunctions[:, ie+1] / 2 - reverse(eigenfunctions[:, ie+1]) / 2

        Norm_sym1 = norm(phi_sym1)
        Norm_sym2 = norm(phi_sym2)

        if Norm_sym1 > Norm_sym2
            Norm_asy_2 = norm(phi_asy2)
            eigenfunctions[:, ie] = phi_sym1 / Norm_sym1
            eigenfunctions[:, ie+1] = phi_asy2 / Norm_asy_2
        else
            Norm_asy_1 = norm(phi_asy1)
            eigenfunctions[:, ie] = phi_sym2 / Norm_sym2
            eigenfunctions[:, ie+1] = phi_asy1 / Norm_asy_1
        end
    end

    return eigenfunctions
end











