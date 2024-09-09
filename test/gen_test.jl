#Test (par affichage) pour la sortie de MSolve
# (field's charasteristic, number of equations, number of variables)
# (7, 2, 2)
function test_msolve_gensys_1()
    d = 7
    neq = 2
    nvar = 2
    eqs, sols, s = GenMSolve(d, nvar, neq)
    display("Equations")
    display(eqs)
    display("MSolve output")
    display(sols)
    display("Calculated solutions")
    display(s)
end

test_msolve_gensys_1()
