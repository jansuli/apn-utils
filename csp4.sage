from tqdm import tqdm
import pickle
from numpy.random import randint
from numpy import array_split
from constraint import *


k=5
n = 2^k-1
m = binomial(n,3) + binomial(n,4)
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
	
basis = [y^i for i in range(2*k)]

def kVecTok2field( vec, offset=0):
	elem = G(0)
	for i in range(len(vec)):
		if vec[i] != 0:
			elem += basis[i+offset]
	#print("Transformed vector %s into field element %s (with offset %d)."%(str(vec), str(elem), offset))		
	return elem

sub = []
for i in range(n):
	sub.append(kVecTok2field(vector(w^i), k))
	

xT = []		# top of columns
xB = []
for i in range(n):
	xT.append(kVecTok2field(vector(w^i), 0))
	xB.append(kVecTok2field(vector(f(w^i)),k))
		
print("Building row transformation matrix 'T'")
		
A = matrix(G, 0, n)
# Generate A
print("Building matrix A. May take some time.")
for ind in tqdm(combInd):
	newRow = matrix(G, 1, n)
	newRow[:,ind] = 1
	A = A.stack(newRow)

x = vector(G, xT) + vector(G, xB)
b= A*x

def getFunc(comp, bivariate=True):
	if bivariate:
		def func(x,y):
			return x+y != comp
	else:
		def func(x):
			return x != comp
	return func
	

def check2Columns(ind1, ind2):
	variables = [ind1, ind2]
	p = Problem()
	for var in variables:
		p.addVariable(var, sub)

	reducedMatrix = A[:, [ind1,ind2]]
	count = 0
	for row in reducedMatrix.rows():
		#print row
		if row == 0 and b[count] == 0:
			print("In position %d impossible.")
			return False
		elif row == 0 and b[count] != 0:
			#print("lucky, no restriction")
			count += 1
		else:
			z = row.nonzero_positions()
			if len(z) == 2:
				fx = getFunc(b[count])
				p.addConstraint(fx, [ind1, ind2])
			else:
				pos = z[0]
				fx = getFunc(b[count], False)
				p.addConstraint(fx, [variables[pos]])
			count += 1
	return p

testIndices = Combinations(n,2).list()
for indexPair in testIndices:
	p = check2Columns(indexPair[0], indexPair[1])
	it = p.getSolutionIter()
	while True:
		try:
			sol = it.next()
			print("Found a soultion: %s"%str(sol))
			with open("cspk=%d_2Cols"%k, "w") as f:
				pickle.dump(sol, f)
		except StopIteration:
			print("Trying different columns")
			break
			
	

		
