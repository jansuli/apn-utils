from tqdm import tqdm
import pickle
from numpy.random import randint
from random import sample
from numpy import array_split
from constraint import *
from os import path
from time import time

k=6
n = 2^k-1
m = binomial(n,3) + binomial(n,4)
print("k = %d, n = %d, m = %d."%(k,n,m))

print("Setting up fields and VectorSpace")

if path.exists("MatrixA_k=%d.sageData"%k):
	with open("MatrixA_k=%d.sageData"%k, "r") as f:
		A = pickle.load(f)
		G = A.base_ring()
		y = G.primitive_element()
else:
	G.<y> = GF(2^(2*k), 'y')
	A = matrix(G, 0, n)

	# Generate A
	print("Generating indices...")
	#comb2 = Combinations(n,2).list()
	comb3 = Combinations(n,3).list()
	comb4 = Combinations(n,4).list()
	combInd = comb3 + comb4
	
	print("Building matrix A. May take some time.")
	for ind in tqdm(combInd):
		newRow = matrix(G, 1, n)
		newRow[:,ind] = 1
		A = A.stack(newRow)
	with open("MatrixA_k=%d.sageData"%k, "w") as f:
		pickle.dump(A, f)

K.<w> = GF(2^k, 'w', repr='log')
R.<z> = PolynomialRing(K,'z')
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

def k2fieldTokField( elem, offset=0):
	res = K(0)
	elemVec = vector(elem)
	for i in range(0,k):
		if elemVec[i + offset] != 0: res += w^i
	return res
	
sub = []
for i in range(n):
	sub.append(kVecTok2field(vector(w^i), k))	

xT = []		# top of columns
xB = []
for i in range(n):
	xT.append(kVecTok2field(vector(w^i), 0))
	xB.append(kVecTok2field(vector(f(w^i)),k))
		
print("First, set of columns...")
Acols = A.columns()
print("Calculating inhomogenity.")

xi = vector(G, xT) + vector(G, xB)
b = vector(G, m)
for i in tqdm(range(n)):
	b = b + xi[i]*Acols[i]

def getFunc(comp, multivariate=True):
	if multivariate:
		def func(*args):
			return sum(args) != comp
	else:
		def func(x):
			return x != comp
	return func
	

def check2Columns(listOfColIndices):
	variables = listOfColIndices[:]
	p = Problem()
	for var in variables:
		p.addVariable(var, sub)

	reducedMatrix = A[:, listOfColIndices]
	count = 0
	print("Adding constraints...")
	constraints = 0
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
			if len(z) == 1:
				pos = z[0]
				fx = getFunc(b[count], False)
				p.addConstraint(fx,[variables[pos]])
			else:
				fx = getFunc(b[count])
				p.addConstraint(fx, variables)
			count += 1
			constraints += 1
	print("Added %d constraints."%constraints)
	return p
	
def updatedXVec(sol):
	XB = xi[:]
	for ind in sol:
		XB[ind] = XB[ind] + sol[ind]
	return XB
	
def polynomialFromVec(apnVec):
	p = R(0)
	image = set()
	for i in range(n):
		res = k2fieldTokField(apnVec[i] + xT[i], k)
		image.add(res)
		p = p + (res*(K(1) + (z+ w^i)^n))
	if len(image) == n:
		print("Found a PERMUTATION")
		with open("PERMUTATIONk=%d_%.0f"%(k, time()), "w") as f:
			f.write(str(p))
	print("It's image has size %d"%len(image))
	return p

print("Looking for solutions...")
testIndices = [sorted(sample(range(n),min(floor(n/2), 15))) for i in range(20)]
sols = 0
solutions = []
for indexPair in testIndices:
	print("testing columns %s"%str(indexPair))
	p = check2Columns(indexPair)
	it = p.getSolutionIter()
	while True:
		try:
			sol = it.next()
			sols += 1
			print("Found a soultion (%d): %s"%(sols,str(sol)))
			apn = updatedXVec(sol)
			pol = polynomialFromVec(apn)
			print("Corresponds to function %s"%str(pol))
			solutions.append(sol)
			if sols < 25:
				with open("cspAPN_k=%d_%d"%(k, sols), "w") as f:
					f.write(str(pol))
		except StopIteration:
			print("Trying different columns")
			break
			
	

		
