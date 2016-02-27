import numpy as np
import numpy.linalg as linalg

def tensorp(A, B):
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
        fu = np.array([[0,1,0,0], [0,0,0,0], [0,0,0,1], [0,0,0,0]], dtype=float)
        fd = np.array([[0,0,1,0], [0,0,0,1], [0,0,0,0], [0,0,0,0]], dtype=float)
        Sc = np.array([[0,1], [0,0]], dtype=float)
        Sd = Sc.T
        fut = fu.T
        fdt = fd.T
        I4 = np.dot(fu, fut) + np.dot(fut, fu)
        I2 = np.dot(Sc, Sd) + np.dot(Sd, Sc)
        Sz = (np.dot(Sc, Sd) - np.dot(Sd, Sc))
        Sx = (Sc + Sd)
        Sy = (Sd - Sc)
        Lx = (np.dot(fut, fd) + np.dot(fdt, fu))
        Ly = (np.dot(fdt, fu) - np.dot(fut, fd))
        Lz = (np.dot(fut, fu) - np.dot(fdt, fd))
        H = - J*(tensorp(Sx, Lx) - tensorp(Sy, Ly) + tensorp(Sz, Lz))
        L2 = Lx**2 - Ly**2 + Lz**2
        S2 = Sx**2 - Sy**2 + Sz**2
        self.I4 = I4
        self.I2 = I2
        self.H =  H
        self.fu_base = tensorp(I2, fu)
        self.fd_base = tensorp(I2, fd)
        self.L2_base = tensorp(I2, L2) + tensorp(S2, I4)
        self.Lz_base = tensorp(I2, Lz) + tensorp(Sz, I4)
        self.lamb = float(lamb)
        self.base = tensorp(I2, I4)
        self.fu = fu
        self.fd = fd
        self.L2 = L2
        self.Lz = Lz
    def __repr__(self):
        return "O GRN  %s. %s"%(self.lamb, self.L2)
    def iteration(self):
        dim = (self.H).shape[0]
        IH = np.identity(dim)
        J = tensorp(self.fu_base, (self.fu).T) + tensorp(self.fd_base, (self.fd).T)
        Jt = J.T
        self.H = tensorp(self.H, self.I4)
        self.H = (self.lamb)*(self.H) + J + Jt
        self.fu_base = tensorp(IH, self.fu)
        self.fd_base =  tensorp(IH, self.fd)
    def diagonalization(self, a):
        if a=="H":
            self.val, self.U = linalg.eigh(self.H)
    def org(self, a):
        if a=="H":
            dim = (self.H).shape[0]
            L = np.zeros((dim))
            for i in range(dim):
                L[i] = i
            for i in range(dim):
                for j in range(dim):
                    if i>j:
                        if self.val[i]<self.val[j]:
                            A = self.val[i]
                            self.val[i] = self.val[j]
                            self.val[j] = A
                            l = L[i]
                            L[i] = L[j]
                            L[j] = l
            I = np.zeros((dim,dim))
            for i in range(dim):
                j = L[i]
                I[j,i] = 1
            self.U = np.dot(self.U, I)
            UI=(self.U).T
            self.H = np.zeros((dim,dim))
            for i in range(dim):
                (self.H)[i,i] = self.val[i] - self.val[0]
            self.fu_base = np.dot(np.dot(UI, (self.fu_base)),self.U)
            self.fd_base = np.dot(np.dot(UI, (self.fd_base)),self.U)
    def cutoff(self, s):
        self.H = (self.H)[0:s,0:s]
        self.fu_base = (self.fu_base)[0:s,0:s]
        self.fd_base = (self.fd_base)[0:s,0:s]
