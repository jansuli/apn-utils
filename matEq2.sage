from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
from tqdm import tqdm
import pickle
import os
from numpy.random import randint


k= 6
n = 2^k-1
m = binomial(n, 2) + binomial(n,3) + binomial(n,4)
print("k = %d, n = %d, m = %d."%(k,n,m))

print("Generating indices...")
combInd = Combinations(n, 4).list() + Combinations(n,3).list() + Combinations(n,2).list()

print("Setting up fields and VectorSpace")
G.<y> = GF(2^(2*k), 'y')
K.<w> = GF(2^k, 'w')
V = VectorSpace(G, n)
VBasis = V.basis()
	
def f(x):
	return x^3
	
basis = [y^i for i in range(2*k)]
def kVecTok2field( vec, offset=0):
	elem = G(0)
	for i in range(len(vec)):
		if vec[i] != 0:
			elem += basis[i+offset]
	print("Transformed vector %s into field element %s (with offset %d)."%(str(vec), str(elem), offset))		
	return elem
	
sub = []
for i in range(n):
	sub.append(kVecTok2field(vector(w^i), k))
	
A = matrix(G, 0, n)
# Generate A
columnOneIndices = dict()

for i in range(n):
	columnOneIndices[i] = set()

print("Building matrix A. May take some time.")
rowCount = 0
for ind in tqdm(combInd):
	newRow = matrix(G, 1, n)
	newRow[:,ind] = 1
	A = A.stack(newRow)
	
	# update dict for later on
	for colIndex in ind:
		columnOneIndices[colIndex].add(rowCount)
	rowCount += 1

print("Testing with gold function:")
xT = []		# top of columns
xB = []		# bottom of columns
for i in range(n):
	xT.append(kVecTok2field(vector(w^i), 0))
	xB.append(kVecTok2field(vector(f(w^i)),k))

xT = vector(G, xT)
xB = vector(G, xB)
x = xT + xB

print 0 in A*x

#b = A*x

b = A*xT
print("Generating constraints")
avoid = set()
for i in range(m):
	row = A[i, :]
	comp = matrix(G,[b[i]])
	rowAug = row.augment(comp)
	if row.rank() == rowAug.rank():
		print("Really a constraint...")
		sol = row.solve_right(comp).columns()[0]
		kernel = row.right_kernel(comp)
		affSpace = AffineSubspace(sol, kernel)
		
		avoid.add(affSpace)
avoid = list(avoid)

def checkConstraints(vec):
	for aff in avoid:
		if vec in aff:
			return False
	else:
		return True

solutions = []
lookedAt = [xB]
# just vary the solution randomly and see what happens
for var in lookedAt:
	print("Checking %s."%str(var))
	for i in range(n):
		for elem in sub:
			test = var[:]
			test[i] = elem
			lookedAt.append(test)
			#print("Testing %s."%str(test))
			if checkConstraints(test) == True:
				print ("Success!!!! Found a solution. We now have found %d apns."%len(solutions))
				solutions.append(test)
			
