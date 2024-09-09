using Mcrypt
using Nemo
using Test

@testset "Mcrypt.jl" begin
    #include("gen_test.jl") #test affichage
    S, x = PolynomialRing(GF(7), ["x", "y"])
    
    for i = 0:6 # tests Issol
        @test Issol([x[1] + i], [7 - i, 0])
        @test Issol([x[2] + x[1] + i], [i, 7 - 2 * i])
    end
    @test Issol([x[1] + 1, x[2] + 2, x[1] * x[2] + x[2] + x[1] + 1], [6, 5])
    
    for i = 1:5 # tests GenSysSol
        eqs, s = GenSysSol(GF(7), ["v", "w", "x", "y", "z"], 5)
        @test Issol(eqs, s) == true
    end
    
    for i = 1:5 # tests GenRand
        eqs, s = GenRand(7, 5, 5, 2)
        @test Issol(eqs, s) == true
    end
    
    for i = 1:5 # tests GenMsolve
        eqs, sols, s = GenMsolve(7, 5, 5, 2)
        @test Issol(eqs, s) == true
        @test Issol(sols, s) == true
    end
    
    for i = 1:5 # tests GenSysPower
        eqs, s = GenSysPower(5, i, 7)
        @test Issol(eqs, s) == true
    end
    
    for i = 1:5 # tests GenSysPowerMsolve
        eqs, sols, s = GenPowerMsolve(5, i, 7)
        @test Issol(eqs, s) == true
        @test Issol(sols, s) == true
    end
    
    # tests PolyToMat
    eqs, s = GenRand(7, 5, 5, 0)
    K = parent(eqs[1])
    ma = [K(1); gens(K)]
    for i in 1:5
        mat = PolyToMat(eqs[i])
        @test transpose(ma) * PolyToMat(eqs[i]) * ma == eqs[i]
    end

    #test GenInstMinRank
    M, lamb = GenInstMinRank(7, 3, 5, 5, 3)
    MSpace = MatrixSpace(GF(7), 5, 5)
    MSum = sum(Int(lamb[i].data) * M[i+1] for i in 1:3) - M[1]
    @test rank(MSpace(MSum)) == 3
end
