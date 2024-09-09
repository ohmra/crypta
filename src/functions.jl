using Nemo
using Random
using AlgebraicSolving

#Generates a vector of quadratic equations using a optimized sum
#Parameters : field, array of variables, number of equations
function GenRandMqSum(K, var, m)
    S, x = PolynomialRing(K, var)
    n = length(x)
    eqs = Array{gfp_mpoly}(undef, m)
    for s in 1:m
        f = sum([rand(K) * x[i] * x[j] for i in 1:n for j in i:n])
        f += sum([rand(K) * x[i] for i in 1:n])
        f += rand(K)
        eqs[s] = f
    end
    eqs
end

#Generates a vector of quadratic equations using matrix product
#Parameters : field, array of variables, number of equations
function GenRandMqMat(K, var, m)
    S, x = PolynomialRing(K, var)
    n = length(x)
    eqs = Array{gfp_mpoly}(undef, m)
    ma = [S(1); x]
    for s in 1:m
        M = rand(K, n+1, n+1) #generate a random matrix
        M = (M + transpose(M))      #symetrize the matrix
        f = transpose(ma) * M * ma
        eqs[s] = f
    end
    eqs
end

#Generates a vector of quadratic equations with its vector solution
#Parameters : field, array of variables, number of equations
function GenSysSol(K, var, m)
  S, x = PolynomialRing(K, var)
  n = length(x)
  s = [rand(K) for i in 1:n]
  eqs = Array{gfp_mpoly}(undef, m)
  p = GenRandMqMat(K, var, m)
  for i in 1:m
    eqs[i] = p[i] - evaluate(p[i], s)
  end
  eqs, s
end

#Intermediary function for easier parameterization
#Parameters : field's characteristic, number of variables
function ParamFieldVar(q, n)
    K = GF(q)
    var = [string("x[", i, "]") for i in 1:n]
    K, var
end

#Generates a vector of quadratic equations using sums / matrices, with / without solution
#Parameters : field's characteristic, number of variables, number of equations, choice of method
function GenRand(q, n, m, method = 1)
    K, var = ParamFieldVar(q, n)
    if method == 0
        eqs = GenRandMqSum(K, var, m)
        s = []
    else
        if method == 1
            eqs = GenRandMqMat(K, var, m)
            s = []
        else
            eqs, s = GenSysSol(K, var, m)
        end
    end
    eqs, s
end

#Generates a system of equations and creates its Groebner basis with MSolve
#Parameters : field's characteristic, number of variables, number of equations
function GenMsolve(q, n, m, method = 2)
    eqs, s = GenRand(q, n, m, method)
    I = AlgebraicSolving.Ideal(eqs)
    sols = groebner_basis(I)
    eqs, sols, s
end

#Generates an instance of minimal rank
#Parameters : field, nb_instance, rows, cols, rank
function GenInstMinRank(q, k, m, n, R)
    K = GF(q)
    lambda = [rand(K) for i in 1:k]
    M = [rand(K, m, n) for i in 1:k]
    MSpace = MatrixSpace(K, m, n)
    A = Matrix(randmat_with_rank(MSpace, R))
    M0 = sum([Int(lambda[i].data) * M[i] for i in 1:k]) - A
    pushfirst!(M, M0) #M_0 = sum(M) - A
    return M, lambda
end


#Generate a vector of equations with v * M = 0
#where M = sum(lambda_i * M_i), M_i generated randomly and v a vector of v_i variables
#Parameters
#q : characteristic of the field, n : (number_of_variables = number_of_equations = number_of_rows = number_of_columns)
function MinRank(q, n)
    lambda_i = [string("lamb[", i, "]") for i in 1:n]
    v_i = [string("v[", i, "]") for i in 1:n]
    K = GF(q)
    S, vars = PolynomialRing(K, [lambda_i; v_i])
    lambda = vars[1:n]
    v = transpose(vars[n+1:n*2])
    M = sum([Matrix(diagonal_matrix(lambda[i], n)) * rand(K, n, n) for i in 1:n])
    eqs = transpose(v * M)
    return eqs
end

#converts a quadratic polynomial to its matrix
#parameter : p, the polynomial of type gfp_mpoly
function PolyToMat(p::gfp_mpoly)
    K = parent(p)
    n = nvars(K)
    vars = [K(1); gens(K)]
    M = zeros(K, n+1, n+1)
    for i in 1:n+1
        for j in i:n+1
            if i==j
                M[i, j] = K(coeff(p, vars[i]*vars[j]))
            else
                M[i, j] = K(coeff(p, vars[i]*vars[j]) // 2)
                M[j, i] = M[i, j]
            end
        end
    end
    return M
end

#Test if vect is a solution of the system of equations
function Issol(syst, vect)
    for eq in syst
        if evaluate(eq, vect) != 0
            return false
        end
    end
    return true
end

#Generate a vector of equations with a solution or not
#n: number_of_variables = number_of_equations, d = degree, q = characteristic of a field
function GenSysPower(n, d, q, solution = true)
    K, var = ParamFieldVar(q, n)
    S, x = PolynomialRing(K, var)
    eqs = Array{gfp_mpoly}(undef, n)
    s = []

    for i in 1:n
        f = 1
        for j in 1:d
            phi_i = sum([rand(K) * x[k] for k in 1:n])
            phi_i += rand(K)
            
            f *= phi_i
        end
        eqs[i] = f
    end

    if solution
        s = [rand(K) for i in 1:n]
        for i in 1:n
            eqs[i] = eqs[i] - evaluate(eqs[i], s)
        end
    end
    eqs, s
end

#Generate a vector of equations with a solution or not and compute its Groebner basis
#n: number_of_variables = number_of_equations, d = degree, q = characteristic of a field
function GenPowerMsolve(n, d, q, solution = true)
    s = []
    if solution
        eqs, s = GenSysPower(n, d, q)
    else
        eqs = GenSysPower(n, d, q, false)
    end
    I = AlgebraicSolving.Ideal(eqs)
    sols = groebner_basis(I)
    eqs, sols, s
end
