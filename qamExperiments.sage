from constraint import *
import multiprocessing as mp
from numpy import array_split
from time import sleep

n = 8
k = 2 # change 2 columns
K.<w> = GF(2^n, 'w' ,repr="log")
Kset = set(K.list())
R.<x> = PolynomialRing(K, 'x')

KBasis = [w^i for i in range(n)]
DualBasis = K.dual_basis()

canonicalBasis = {
	n-1 : VectorSpace(GF(2), n-1).basis(),
	n-2 : VectorSpace(GF(2), n-2).basis(),
	}

combinations = dict()

for i in range(n+1):
	combinations[i] = Combinations(i).list()

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
		#print N
		for ind in combinations[N]:
			#print ind
			ret.add( K(sum([row[i] for i in ind])) )
		return ret	
def f(x):
	return x^3
	
Cf = matrix(K,n,n)
Cf [0,1] = Cf[1,0] = 1

H = M.transpose()*Cf*M
print H.str() + "\n"

A = H[:n-1, :n-1]

KSet = set(K.list())

########## Treat as CSP ########

def Sf(indexVec, mat):
	return KSet.difference( rowSpan(indexVec*mat))

def getDomains(mat):
	domains = dict()
	N = mat.nrows()
	for i in range(N):
		domains[i] = list(Sf(canonicalBasis[N][i], mat))
	return domains
		
def getConstraintFunc(avoid):
	def SumNotInSet(*args):
		return sum(args) not in avoid
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
		constraintList.append((constraintFunction,ind))

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

newCol = matrix(K, n-1, 1)
S = dict()
for ind in combinations[n-1]:
	if len(ind) > 0:
		S[tuple(ind)] = Kset.difference(rowSpan( sum(A[ind, :].rows())))
		
def getSetFn(oldFn, xi, pos):
	def newSetFn(indexTuple):
		return oldFn(indexTuple).intersection(set([xi + v for v in oldFn(tuple([pos]+list(indexTuple)))]))
	return newSetFn	
		
def nextComponent(columnTilNow, setFn):
	newCol = copy(columnTilNow)
	pos = n - (1+newCol.list().count(0))
	domain = setFn((pos,))
	if len(domain) > 0:
		print("At %d domain has length %d."%(pos, len(domain)))
		for xi in domain:
			newCol[pos] = xi			
			if pos < A.nrows()-1:
				newSetFn = getSetFn(setFn, xi, pos)
				next = nextComponent(newCol,newSetFn)
				if next != None: return next
			else:
				return newCol
	
def setFn( indexTuple):
	return S[indexTuple]
	
col = nextComponent(newCol, setFn)
