# NRG-Didactic

Modules that runs on Python with Numpy that define classes called NRG. This classes has all the attributes and methods necessary to make a Numerical Renormalization Group (NRG). There is two modules: 'Kondo.py' and 'Coupled-Spins.py'. The first is the NRG applied to Kondo Model and the second is an toy model of spins coupled with the neighborhood by an exponential decreasing function "$J\lambda^ n$"

#Initiating the class

The module 'Kondo.py' is initiated by the parameters
- Lambda parameter of logarithm discretization; lambda
- Coupling of the spin impurity with the spins of the electrons; J

The next module called 'coupled-spins.py' is initiated by the parameters
- Base of the exponential decreasing function; lambda
- Front constant of the exponential decreasing function; J

#Methods ofthe Class 

There are three main methods: diagonalization of the hamiltonian '.diag()', iteration of the hamiltonian '.ite()' and the .cutoff(S).
- .diag(): This method change the basis of all matrices that are important for the '.ite()', diagonalizing the hamiltonian
- .ite(): This method apply the numerical renormalization group transformation.
- .cutoff(S): This method cutoff the higer excited states, holding only the first S states excitate. (S is int)

There are lot of attributes that works internilly in the class. The important attribute are the hamiltonian '.H'(a two dimensional numpy array) and the energies '.val'(a one dimensional numpy array).

#Example of an NRG process for the 'Kondo.py' module:

-nrg = NRG(1.41421356, -1)                       #Initialization of the class#
-E = np.zeros((4000,4000))                       #An 2 dimensional array to record the eigenvalues of any iteration#
-for i in range(200):
-    print ("Diagonalization %s"%i)
-    nrg.diag()                                  #Diagonalization. Here the nrg.H is not only diagonalized but rescaled too#
-    if nrg.val.shape[0] > 1024:
-          nrg.cutoff(1024)                      #Cut off the higest states kepting only with 1024 first states#
-    print ("Iteration %s"%i)
-    nrg.ite()                                   #Construction of the new matrix#
-    print("Memory %s"%i)
-    for j in range(nrg.val.shape[0]):
-	      E[i,j] = nrg.val[j] - nrg.val[0]        #Recording the eigenvalues of the 'i' iteration rescaled#
-    print ("-------------")



