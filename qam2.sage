from simpleai.search import CspProblem
from simpleai.search import backtrack, min_conflicts, MOST_CONSTRAINED_VARIABLE, LEAST_CONSTRAINING_VALUE, convert_to_binary
import multiprocessing as mp
from numpy import array_split
from time import sleep

n = 7
k = 2 # change 2 columns
K.<w> = GF(2^n, 'w' ,repr="log")
R.<x> = PolynomialRing(K, 'x')

KBasis = [w^i for i in range(n)]
DualBasis = K.dual_basis()

canonicalBasis = {
	n-1 : VectorSpace(GF(2), n-1).basis(),
	n-2 : VectorSpace(GF(2), n-2).basis(),
	}

combinations = {
	n-2 : Combinations(n-2).list(),
	n-1 : Combinations(n-1).list(),
	n : Combinations(n).list(),	# nontrivial row combinations for qam checking
	}

M = matrix(K, n,n)
for i in range(n):
	for j in range(n):
		M[i,j] = KBasis[j]^(2^i)
		
MB = matrix(K,n,n)
for i in range(n):
	for j in range(n):
		MB[i,j] = DualBasis[j]^(2^i)
		
def checkRowRank(row):
	testMatrix = matrix( GF(2), n,0)
	for i in range(row.ncols()):
		testMatrix = testMatrix.augment(matrix(vector(row[0,i])).transpose())
	#print testMatrix.str()
	return testMatrix.rank()

def checkQAM(mat):
	N = mat.ncols()
	for ind in combinations[N][1:]:		
		testGenerator = matrix(sum(mat[ind,:]))
		#print testGenerator
		if checkRowRank(testGenerator) != N-1:
			print("no qam caus of sum of %s"%str(ind))
			return False
	else:
		print ("QAM")
		return True

def rowSpan(row, matrix = False):
	ret = set()
	if matrix:
		N = row.ncols()
		print "Matrix"
		for ind in combinations[N]:
			ret.add( sum( row[0, ind].list() ))
		return ret
	else:
		N = len(row)
		for ind in combinations[N]:
			ret.add( sum([row[i] for i in ind]) )
		return ret	
	
def f(x):
	return x^3
	
Cf = matrix(K,n,n)
Cf [0,1] = Cf[1,0] = 1

H = M.transpose()*Cf*M
print H.str() + "\n"

A = H[:n-2, :n-2]

KSet = set(K.list())

########## Treat as CSP ########

def Sf(indexVec, mat):
	return KSet.difference( rowSpan(indexVec*mat))

def getDomains(mat):
	domains = dict()
	N = mat.nrows()
	print N
	for i in range(N):
		print i
		domains[i] = list(Sf(canonicalBasis[N][i], mat))
	return domains
		
def getConstraintFunc(avoid):
	def SumNotInSet(variables, values):
		#print values
		return sum(values) not in avoid
	return SumNotInSet 
	
def	getConstraints(mat):
	'''Takes submatrix and generates constraints such that no sum of variables might lay in the row-span of the corresponding sum of rows.'''
	
	N = mat.ncols()
	sumIndices = [ind for ind in combinations[N] if len(ind) > 1]

	constraintList = []
	for ind in sumIndices:
		indexVec = matrix(GF(2), 1, N)
		indexVec[0, ind] = 1
		avoidSet = list(rowSpan((indexVec*mat).list()))
		#print indexVec, avoidSet
		constraintFunction = getConstraintFunc(avoidSet)
		constraintList.append((tuple(ind), constraintFunction))

	return constraintList

def applySolution(mat, solDict):
	N = max(solDict.keys()) + 1
	#print("MAX IS:%d"%N)
	for i in solDict:
		mat[N, i] = mat[i, N] = solDict[i]
		
def solutionPolynomial(qam):
	C = MB*qam*MB.transpose()
	p = R(0)
	for i in range(1,qam.ncols()):
		for j in range(i):
			p+= C[i,j]*x^(2^i + 2^j)
	return p
	
def iterate(it, oldMat, k, polList):
	while True:
		try:
			sol = it.next()
			Hmod = copy(oldMat)
			applySolution(Hmod,sol)
			print ("Found a column such that %d×%d submatrix is qam:\n%s\n"%(n-k+1, n-k+1,Hmod.str()))
			pol = solutionPolynomial(Hmod)
			#print Hmod[:n-k+1, :n-k+1]
			#checkQAM(Hmod[:n-k+1, :n-k+1])
			print( "The corresponding polynomial is %s."%str(pol))
			
			if k == 2:
				print("Now looking for alternative next columns.")
				cspWorker(Hmod, k-1, polList)
			elif k == 1:
				polList.append(pol)
		except StopIteration:
			print ("No more solutions.")
			return
		
def cspWorker(mat, k, polList, domains = None):
	A = mat[:n-k, :n-k]
	print("Setting up problem.")
	print("Generating domains and adding variables.")
	if not domains:
		print ("No domains given...")
		domains = getDomains(A)
	variables = domains.keys()
	print("Now adding constraints.")
	constraints = getConstraints(A)
	
	print("Converting to binary constraint problem.")
	variables, domains, constraints = convert_to_binary(variables, domains, constraints)
	
	p = CspProblem(variables, domains, constraints)
	print("%d vars and %d constraints added."%(len(domains), len(constraints)))
	print("Starting iteration.")
	result = backtrack(p,variable_heuristic=MOST_CONSTRAINED_VARIABLE,value_heuristic=LEAST_CONSTRAINING_VALUE)
	return result
	
# Starting multicore-search
pols = []
sol = cspWorker(H, k, pols)
print("Found solution for reduced problem.")
applySolution(H, sol)
print("Now H is \n%s"%H.str())
sol2 = cspWorker(H, k-1, pols)
applySolution(H, sol2)
