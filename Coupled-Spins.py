import numpy as np
import numpy.linalg as linalg

def tensorp(A, B):
    cdef int i, j, k, l, dim, dimA, dimB, dgim, ni, nj
    dimA=A.shape[0]
    dimB=B.shape[0]
    gdim=dimA*dimB
    G=np.zeros((gdim, gdim))
    for l in range(dimB):
        for k in range(dimB):
            for i in range(dimA):
                for j in range(dimA):
                    ni=(dimB*i)+k
                    nj=(dimB*j)+l
                    G[ni,nj]=B[k,l]*A[i,j]
    return G

class NRG(object):
    def __init__(self, lamb, J):
        Sc = np.array([[0,1], [0,0]], dtype=float)
        Sd = Sc.T
        I2 = np.dot(Sc, Sd) + np.dot(Sd, Sc)
        Sz = (np.dot(Sc, Sd) - np.dot(Sd, Sc))
        Sx = (Sc + Sd)
        Sy = (Sd - Sc)
        H0 = - (J)*(tensorp(Sx, Sx) - tensorp(Sy, Sy) + tensorp(Sz, Sz))
        self.I2 = I2
        self.H =  H0
        self.Sc_base = tensorp(I2, Sc)
        self.Sd_base = tensorp(I2, Sd)
        self.lamb = float(lamb)
        self.J = J
        self.Sx = Sx
        self.Sy = Sy
        self.Sz = Sz
        self.Sc = Sc
    def __repr__(self):
        return "O GRN  %s. %s"%(self.lamb, self.J)
    def ite(self):
        dim = (self.H).shape[0]
        IH = np.identity(dim)
        Sx_base = self.Sc_base + self.Sd_base
        Sy_base = self.Sd_base - self.Sc_base
        Sz_base = np.dot(self.Sc_base, self.Sd_base) - np.dot(self.Sd_base, self.Sc_base)
        V = - (self.J)*(tensorp(Sx_base, self.Sx) - tensorp(Sy_base, self.Sy) + tensorp(Sz_base, self.Sz))
        self.H = tensorp(self.H, self.I2)
        self.H = (self.lamb)*(self.H) + V
        self.Sc_base = tensorp(IH, self.Sc)
    def diag(self):
        self.val, self.U = linalg.eigh(self.H)
        dim = self.U.shape[0]
        UI=(self.U).T
        self.H = np.zeros((dim,dim))
        for i in range(dim):
            (self.H)[i,i] = self.val[i] - self.val[0]
        self.Sc_base = np.dot(np.dot(UI, (self.Sc_base)),self.U)
        self.Sd_base = (self.Sc_base).T
    def cutoff(self, s):
        self.H = (self.H)[0:s,0:s]
        self.Sc_base = (self.Sc_base)[0:s,0:s]
        self.Sd_base = (self.Sc_base).T