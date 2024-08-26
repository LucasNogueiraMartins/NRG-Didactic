# NRG-Didactic

These modules run on Python using the Numpy library. It defines a class called NRG. This class has all the attributes and methods necessary for implementing the Numerical Renormalization Group (NRG). There are two modules: 'Kondo.py' and 'Coupled-Spins.py'. The first is the NRG applied to Kondo Model and the second is a toy model of spins with neighborhood coupling and an exponential decreasing function "$J\lambda^ n$"

#Initiating the class

The module 'Kondo.py' is initiated by the parameters
- Parameter that sets the logarithm discretization; lambda
- Coupling between the spin impurity and the spins of the electrons; J

The next module called 'coupled-spins.py' is initiated by the parameters
- Base of the exponential decreasing function; lambda
- pre-factor of the exponential decreasing function; J

#Methods of the Class 

There are three main methods: diagonalization of the Hamiltonian '.diag()', iteration of the Hamiltonian '.ite()', and the .cutoff(S).
- .diag(): This method changes the basis of all matrices that are important for the '.ite()', diagonalizing the Hamiltonian
- .ite(): This method applies the numerical renormalization group transformation.
- .cutoff(S): This method cuts off the higher excited states, holding only the first S states excited. (S is int)

There are a lot of attributes that work internally in the class. The important attributes are the Hamiltonian '.H'(a two-dimensional numpy array) and the energies '.val'(a one-dimensional numpy array).


