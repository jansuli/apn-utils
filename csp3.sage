from tqdm import tqdm
import pickle
from numpy.random import randint
from numpy import array_split
from constraint import *


k=4
n = 2^k-1
m = n + binomial(n,3) + binomial(n,4)
print("k = %d, n = %d, m = %d."%(k,n,m))

print("Generating indices...")
#comb2 = Combinations(n,2).list()
comb3 = Combinations(n,3).list()
comb4 = Combinations(n,4).list()

combInd = comb3 + comb4
print("Setting up fields and VectorSpace")
G.<y> = GF(2^(2*k), 'y')
K.<w> = GF(2^k, 'w')
V = VectorSpace(G, m)
VBasis = V.basis()
	
def f(x):
	return x^3
	
print("Building row transformation matrix 'T'")

T = matrix.identity(G,m)
row = n

for ind in tqdm(combInd):
	for addRow in ind:
		T[row] = T[row] + T[addRow]
	row += 1
		
A = matrix.identity(G,n)
# Generate A
print("Building matrix A. May take some time.")
for ind in tqdm(combInd):
	newRow = matrix(G, 1, n)
	newRow[:,ind] = 1
	A = A.stack(newRow)
	
Tsub = T[n:, :]
kern = Tsub.right_kernel()

variables = list(range(n))
p = Problem()
p.addVariables(variables, G.list()[1:])

print("adding constraints")
for baseVec in kern.basis_matrix().transpose().rows():
	affected = [variables[i] for i in baseVec.nonzero_positions()]
	p.addConstraint(lambda *args: sum(args) != 0, affected)
	
def getBVec(sol):
	b = V.linear_combination( (kern.basis()[i], sol.values()[i]) for i in range(n) )
	return b

		
