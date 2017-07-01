from tqdm import tqdm
import pickle
from numpy.random import randint
from numpy import array_split
from constraint import *


k=5
n = 2^k-1
m = binomial(n, 2) + binomial(n,3) + binomial(n,4)
print("k = %d, n = %d, m = %d."%(k,n,m))

print("Generating indices...")
#comb2 = Combinations(n,2).list()
comb3 = Combinations(n,3).list()
comb4 = Combinations(n,4).list()

combInd = comb3 + comb4
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
	#print("Transformed vector %s into field element %s (with offset %d)."%(str(vec), str(elem), offset))		
	return elem
	
sub = [G(0)]
for i in range(n):
	sub.append(kVecTok2field(vector(w^i), k))
	
A = matrix(G, 0, n)
# Generate A

#print("Building matrix A. May take some time.")
#for ind in tqdm(combInd):
	#newRow = matrix(G, 1, n)
	#newRow[:,ind] = 1
	#A = A.stack(newRow)
	
	
#print("Testing with gold function:")
#xT = []		# top of columns
#xB = []		# bottom of columns
#for i in range(n):
	#xT.append(kVecTok2field(vector(w^i), 0))
	#xB.append(kVecTok2field(vector(f(w^i)),k))

#xT = vector(G, xT)
#xB = vector(G, xB)
#x = xT + xB

#print 0 in A*x

##b = A*x

#b = A*xT
p = Problem(MinConflictsSolver())

xT = []		# top of columns
xB = []
for i in range(n):
	xT.append(kVecTok2field(vector(w^i), 0))
	xB.append(kVecTok2field(vector(f(w^i)),k))
	
print("Generating variables for CSP")
variables = []
for i in range(n):
	var = i
	variables.append(var)
	indices = randint(0, len(sub), 25)
	dom = [xB[i] + sub[j] for j in indices]
	p.addVariable(var, dom)

print("Adding constraints...")
for comb in combInd:
	affected = [variables[j] for j in comb]
	compareTop = [xT[j] for j in comb]
	def getFunction(comps):
		def func(*args):
			return sum(comps) + sum(args) != 0		
		return func
	fx = getFunction(compareTop)
	p.addConstraint(FunctionConstraint(fx), affected)
# i=0
#for comb in combInd:
	#affectedVars = [vars[j] for j in comb]
	#def getFunction(comp):
		#def func(*args):
			#return sum(args)!= comp
		#return func
	#fx = getFunction(b[i])
	#p.addConstraint(FunctionConstraint(fx), affectedVars)
	#i += 1
	
def check3Dependence(sol):
	s = []
	sortedSol = [elem[1] for elem in sorted(sol.items())]
	for comp in sortedSol:
		s.append(comp)
	v = [xT[i] + s[i] for i in range(n)]
	for ind in comb3:
		test = [v[i] for i in ind]
		if sum(test) == 0:
			#print("found three lin dep comps in:%s"%str(test))
			return False
	else:
		return True
			
def check4Dependence(sol):
	s = []
	sortedSol = [elem[1] for elem in sorted(sol.items())]
	for comp in sortedSol:
		s.append(comp)
	v = [xT[i] + s[i] for i in range(n)]
	for ind in comb4:
		test = [v[i] for i in ind]
		if sum(test) == 0:
			return False
	else:
		return True
			
print("Done, now looking for solutions")
#it = p.getSolutionIter()
#count = 0
#while True:
	#try:
		#solution = it.next()
		#count += 1
		#print("Saving solution.")
		#print("We currently got %d solutions."%count)
		#with open("csp%dSol%dBottomDict"%(k,count), "w") as f:
			#pickle.dump(solution, f)
	#except StopIteration:
		#break
	

	

		
