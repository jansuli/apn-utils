from simpleai.search import CspProblem
from simpleai.search import backtrack, min_conflicts, MOST_CONSTRAINED_VARIABLE, LEAST_CONSTRAINING_VALUE, convert_to_binary
import multiprocessing as mp
from numpy import array_split
from time import sleep
#from sympy import var as syVar
from sympy.logic.boolalg import to_cnf, conjuncts
from tqdm import tqdm

n = 7
k = 2 if n>8 else 1# change 2 columns
K.<w> = GF(2^n, 'w' )#,repr="log")
R.<x> = PolynomialRing(K, 'x')

KBasis = [w^i for i in range(n)]
DualBasis = K.dual_basis()

combinations = {
	n-3: Combinations(n-3).list(),
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

A = H[:n-2, :n-2]

KSet = set(K.list())

def affineTranslations(rowVector):
	V = rowSpan(rowVector)
	domain = KSet.difference(V)
	translations = []
	build = set()
	for a in domain:
		v = V.pop()
		if a+v not in build:
			V.add(v)
			coset = set( [a+v for v in V] )
			build = build.union(coset)
			translations.append(a)	
		else:
			V.add(v)	
			
	return translations
	
def affineUnion(rowVector, offset = 0):
	print ("affine union")
	translations = affineTranslations(rowVector)
	
	if offset == 0:
		return translations
	else:
		ret = set()
		V = rowSpan( rowVector[offset:] )
		for a in translations:
			ret = ret.union( set([a + v for v in V]))
		return list(ret)

########## Treat as CSP ########

def getDomains(mat):
	domains = dict()
	N = mat.nrows()
	print N
	for i in range(2):
		domains[i] = affineUnion(mat.rows()[i],i)
	for i in range(2,N):
		domains[i] =  list(KSet.difference( rowSpan(mat.rows()[i]) ))
	return domains
		
def getConstraintFunc(avoid):
	def SumNotInSet(variables, values):
		print values
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
		constraintList.append((tuple(ind), avoidSet))

	return constraintList
	
A = H[:n-1, :n-1]

# Build sympy expressions
## not in set constraints
expr = []
for i in range(n-1):
	l = i*n + 1
	rowSpace = rowSpan(A.rows()[i])
	for s in rowSpace:
		clause = "Or( "
		sVec = vector(s)
		k = 0
		for j in sVec:
			if j == 1:
				clause += " ~x%d ,"%(l+k)
			else:
				clause += " x%d ,"%(l+k)
			k += 1
		#print clause[:-1]
		expr.append( clause[:-1] + ")")
				
## in set constraints
#print("\nIn Set\n")
#domains = getDomains(A)
#expr2 = "Or(" 
#for i in domains:
	#l = i*n
	#dom = domains[i]
	#for s in dom:
		#clause = ""
		#sVec = vector(s)
		#k = 0
		#for j in sVec:
			#if j == 1:
				#clause += "x[%d]&"%(l+k)
			#else:
				#clause += "x[%d]&"%(l+k)
			#k += 1
		#print clause
		#expr2 += clause[:-1] + ","
#expr2 = expr2[:-1] + ")"

# sum not in set
constraints = getConstraints(A)
expr3 = []
for ind,avoid in tqdm(constraints, desc="constraints"):
	for s in tqdm(avoid):
		sVec = vector(s)
		k = 0
		clause = ""
		for j in sVec:
			comps = [i*n+k+1 for i in ind]
			comps = ["x%d"%i for i in comps]
			if j == 1:
				clause += "Not(Xor(" + ",".join(comps) + ")) |"
			else:
				clause += "Xor(" + ",".join(comps) + ") |"
			k += 1
		elems = [str(e) for e in conjuncts(to_cnf(clause[:-1]))]
		#tqdm.write(str(elems))
		expr3.append(elems)

cnf = ""
clausesN = 0
for line in expr:
	line = line.replace("Or(", "").replace(")","").replace("x","").replace("~", "-").replace(",","") +" 0\n"
	print(line)
	cnf += line
	clausesN += 1
	
for clauses in expr3:
	print("\nNew list of clauses.")
	for line in clauses:
		line = line.replace("Not(","-").replace("Or(", "").replace("x","").replace(")","").replace(",","") + " 0\n"
		print line
		cnf += line
		clausesN += 1
		
cnf = "p cnf %d %d"%(n*(n-1), clausesN) + "\n" + cnf 

with open("sat%d.cnf"%n ,"w") as f:
	f.write(cnf)
	pass
	
	

		
	
