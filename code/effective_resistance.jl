#=
Code for Sparsification.
Right now, just implements Spielman-Srivastava
=#

using Laplacians
using LinearAlgebra
using MAT

"""
as = sparsify(a; ep=0.5)
Apply Spielman-Srivastava sparsification: sampling by effective resistances.
`ep` should be less than 1.
"""
function resistance(a, pairs; ep=0.0, matrixConcConst=4.0, JLfac=4.0)
    if ep > 1
    @warn "Calling sparsify with ep > 1 can produce a disconnected graph."
    end

    f = approxchol_lap(a,tol=1e-2);
    n = size(a,1)
    k = round(Int, JLfac*log(n))

    @show n

    # number of dims for JL
    U = wtedEdgeVertexMat(a)
    m = size(U,1)
    R = randn(Float64, m,k)
    UR = U'*R;
    V = zeros(n,k)
    for i in 1:k
        V[:,i] = f(UR[:,i])
    end

    """
    rs = zeros(n,n)
    for i in 1:n for j in i:n
        rs[i,j] = (norm(V[i,:]-V[j,:])^2/k)
            # calculate the effective resistance between the ith and jth nodes
    end
    """
    pair_size = size(pairs, 1)
    rs = zeros(pair_size)
    for i in 1:pair_size
       v1 = pairs[i,1] + 1
       v2 = pairs[i,2] + 1
       #@show v1, v2
       rs[i] = (norm(V[v1,:]-V[v2,:])^2/k)
    end
    return rs
end


# There are 3 args: current site index, dataset name, and total number of site
# Read the adjacency matrix and the paris need to calculate the effective resistance
filename1 = string("../data/", ARGS[1], "/tmp/sp_A_n", ARGS[2] , "_", ARGS[3], ".mat")
file = matopen(filename1)
A_sp = read(file, "A_sp")

filename2 = string("../data/", ARGS[1], "/tmp/pairs_n", ARGS[2], "_", ARGS[3], ".mat")
file = matopen(filename2)
pairs = read(file, "pairs")

# calculate the effective resistance for those pairs
time1 = time()
rs = resistance(A_sp, pairs; ep=0.3, matrixConcConst=4.0, JLfac=4.0)
time2 = time()
time_diff = time2-time1
@show time_diff

# save the calculated effective resistance
filename3 = string("../data/", ARGS[1], "/tmp/rs_n", ARGS[2], "_", ARGS[3],".mat")
file = matopen(filename3, "w")
write(file, "rs", rs)
close(file)