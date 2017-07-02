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

#if path.exists("MatrixA_k=%d.sageData"%k):
	#print("Opening")
	#with open("MatrixA_k=%d.sageData"%k, "r") as f:
		#A = pickle.load(f)
		#G = A.base_ring()
		#y = G.primitive_element()
	#print("Done")

G.<y> = GF(2^(2*k), 'y')
combInd = Combinations(n,3).list() + Combinations(n,4).list()
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
	
sub = [G(0)]
for i in range(n):
	sub.append(kVecTok2field(vector(w^i), k))	

xT = []		# top of columns
xB = []
for i in range(n):
	xT.append(kVecTok2field(vector(w^i), 0))
	xB.append(kVecTok2field(vector(f(w^i)),k))
	
xi = vector(G, xT) + vector(G, xB)
		

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
		p.addVariable(var, [xT[var] + s for s in sub if s != xB[var]])

	cols = set(listOfColIndices)
	print("Collecting relevant indices.")
	reducedIndices = [set(ind) for ind in combInd if set(ind).intersection(cols) != set()]
	count = 0
	print("Adding constraints...")
	constraints = 0
	for ind in reducedIndices:
		rhs = G(0)
		for j in ind.difference(cols):
			rhs += xi[j]
		affected = list(cols.intersection(ind))
		fx = getFunc(rhs)
		p.addConstraint(fx, affected)
		constraints += 1		
		
	print("Added %d constraints."%constraints)
	return p
	
def updatedXVec(sol):
	XB = xi[:]
	for ind in sol:
		XB[ind] = sol[ind]
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
			
	

		
