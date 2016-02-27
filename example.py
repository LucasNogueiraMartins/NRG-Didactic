nrg = NRG(2, -0.024) #initiating the class#
E = np.zeros((4000,4000)) #Tabula-rasa#

N = 10   # number of iterations
M = 1024 # number of states kept
for j in range(N):
    print("diagonalization")
    nrg.diag()
    print("memory")
    for i in range((nrg.val).shape[0]):
        E[j, i] = nrg.val[i] - nrg.val[0]
    print("cutoff")
    if nrg.val[0] > M:
        nrg.cutoff(M)
    print("iteration"")
    nrg.ite()
    print"-------------------"