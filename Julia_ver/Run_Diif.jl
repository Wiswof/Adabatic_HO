using Statistics
using LinearAlgebra
using Plots
using ArnoldiMethod, GenericSchur, DoubleFloats 
using DoubleFloats.Quadmath
using Arpack,SparseArrays
setprecision(128)


include("Functions_.jl")


# Constants in atomic units
L = big"5.0" * big"1e2"    # Region size  big"50.0"
hbar = big"1.0"             # Reduced Planck's constant
m = big"1.0"                # Electron mass
L = Float128(L)
hbar =  Float128(hbar)
m=Float128(m)
# Define the size of the region and the number of points for discretizationc
N = Int(1.0e1)   #1.0e6             # Number of points for discretization (5*1e5 is enough for L = 5*1e3)
N_state = 8                 # Number of states from vacuum
dx = L / N                  # Discretization step
x = range(-L / 2, stop=L / 2, length=N)

# Define the potential function parameters
# Vpar = Dict(
#     "Omega" => big"1.5e5" * π,     # Hz
#     "kappa" => big"2.0e6" * π,      # Hz
#     "Delta" => big"2.0e2" * π,      # Hz
#     "g" => big"1875.0" * 2π
# )

Vpar = Dict(
    "Omega" => Float128(big"1.5e5" * π),     # Hz
    "kappa" => Float128(big"2.0e6" * π),      # Hz
    "Delta" => Float128(big"2.0e2" * π),      # Hz
    "g" => Float128(big"1875.0" * 2π)
)
Vpar["nc"]=Vpar["Omega"]*(4*Vpar["Delta"]^2
    +Vpar["kappa"]^2)/(16*Vpar["g"]^2*Vpar["Delta"]);



FV_alpha(alpha) = 0.5 * Vpar["Omega"] * x .^ 2 +
                  alpha * Vpar["nc"] * (Vpar["Delta"]^2 + Vpar["kappa"]^2 / 4) / (Vpar["kappa"] / 2) *
                  (atan.((Vpar["g"] .* x .+ Vpar["Delta"]) / (Vpar["kappa"] / 2)) +
                   atan.((-Vpar["g"] .* x .+ Vpar["Delta"]) / (Vpar["kappa"] / 2)))


plot(x,FV_alpha(1))
WaveThre = 0.01  # truncate value on wavefunction for |psi| < max(|psi|) * WaveThre
Nalpha = Int(100)
Dalpha = Float128(big"0.001")
ial = range(1 - Dalpha, stop=1 + 2 * Dalpha, length=Nalpha)

Dn_nc = mean(diff(ial))


 eigenfunctions, eigenvalues = Eigen_state_value(N, FV_alpha(ial[1] - Dn_nc), hbar, m, dx, N_state)
# Gound_state_ini = eigenfunctions[:, 1]  # N_state copies



# Dphi_thre = 1e-4  # threshold for <Dphi|phi> is between Dphi_thre and Dphi_thre/2
