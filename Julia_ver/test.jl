using ArnoldiMethod, GenericSchur, DoubleFloats ,LinearAlgebra
using DoubleFloats.Quadmath
using Arpack,SparseArrays
T = Float128
#Complex{Double64}

N = 64
64

A = Tridiagonal(rand(T, N - 1), rand(T, N), rand(T, N - 1));
B=sparse(A)
B=eigen(A) 

pschur,hist = partialschur(A; nev=round(Int, 0.1 * N), which=LM());



# 创建一个稀疏矩阵
sparse_matrix = sparse([1, 2, 3], [1, 2, 3], [1.0, 2.0, 3.0])

# 将稀疏矩阵转换为一般矩阵
dense_matrix = Matrix(sparse_matrix)

# 打印结果
println("Sparse Matrix:")
println(sparse_matrix)

println("\nDense Matrix:")
println(dense_matrix)



NN=1000;
A=Float128
A = rand(Float128, 2, 2)#randn(NN, NN)
A = A + A'  # 使矩阵对称
#Float128(randn(2,2))
# 计算特征值和特征向量
v = eigs(A, nev = 1, which=:LM)