#from sympy import var as syVar
from sympy.logic.boolalg import to_cnf, conjuncts
from tqdm import tqdm
import pickle

n = 8
K.<w> = GF(2^n, 'w' ,repr="log")
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

		
def differConstraint(vec):
	'''Takes an n-1 vector of F_{2^n} elements and returns a boolean function as CNF expression in (n-1)*n variables ensuring that the binary vectors differ.'''
	expression = "And( "
	count = 0
	for elem in vec:
		subexpression = "Or( "
		binary = vector(elem)
		print (elem,binary)
		for entry in binary:
			if entry == 1:
				subexpression += "~{%d},"%count
			else:
				subexpression += "{%d},"%count
			count += 1
		expression += subexpression[:-1] + "),"
	expression = expression[:-1] + ")"
	return expression
	
def notInSetConstraint(outSet):
	expression = "And( "
	for elem in outSet:
		subexpression = "Or( "
		vec = vector(elem)
		count = 0
		for entry in vec:
			if entry == 1:
				subexpression += "~{%d},"%count
			else:
				subexpression += "{%d},"%count
			count += 1
		expression += subexpression[:-1] + "),"
	expression = expression[:-1] +")"
	print ("Added %d clauses."%len(outSet))
	return expression
	

def generatePreSat(mat, filename, differ = None):
	'''Takes a quadratic QAM-submatrix 'mat' and dimension/degree of galois field and generates a boolean expression wich has to be further transformed using Tseitin.'''
	N = mat.nrows()
	expressionList = []
	k = mat.base_ring()
	dim = k.degree()
	Kset = set(k)
	
	print("Generating sat problem for submatrix with %d rows over field with degree %d"%(N,dim))
	
	if differ:
		differExpr= differConstraint(differ).format(*['x%d'%i for i in range(1, N*dim+1)])
		print differExpr
		expressionList.append(differExpr)
	
	# Add domain constraints. Domains have size at most 3*2^(N-1), while the set difference has size 2^(N-1).
	# Both constraint types are equivalent and produce same amount of clauses
	for i in range(N):
		if i <= 1:
			aff = affineUnion(mat.rows()[i],i)
			print aff
			rowSpace = Kset.difference(aff)
			#rowSpace = rowSpan(mat.rows()[i])
		else:
			rowSpace = rowSpan(mat.rows()[i])
			print (0 in rowSpace)
		offset = i*dim + 1
		
		formatting = ["x%d"%j for j in range(offset, offset+dim)]
		expr = notInSetConstraint(rowSpace) 
		constraintExpression = expr.format(*formatting)
		expressionList.append(constraintExpression)
		
	expression1 = "And("+ ",".join(expressionList)+")"  # everything up till now is already in CNF and doesn't need to be transformed
	expressionList = []
	print expression1
	sumIndices = [ind for ind in combinations[N] if len(ind) > 1]
	for ind in sumIndices:
		formatting = ["Xor( " for i in range(dim)] # put #dim xor's into constraintFunctions  
		# Iterate over vectors/rows to add:
		for i in ind:
			offset = i*dim + 1
			for j in range(dim):
				formatting[j] += "x%d,"%(offset+j)
		formatting = [formatting[i][:-1] + " )" for i in range(dim)]
		
		rowSum = sum(mat[ind, :].rows()) 
		expr = notInSetConstraint( rowSpan(rowSum) )
		constraintExpression = expr.format(*formatting)
		
		expressionList.append(constraintExpression)
		
	#completeExpression = "And( " + ",".join(expressionList) + " )"
	
	with open(filename, "w") as f:
		pickle.dump((expression1, expressionList),f)
		
	print ("All done and saved under %s."%filename)
	
def applySatSolution(filename, sub):
	'''Reads sat result file and appends a new column to given submatrix.'''
	N = sub.nrows()
	k = sub.base_ring()
	w = k.primitive_element()
	dim = k.degree()
	
	with open(filename, "r") as f:
		for line in f:
			if line.startswith("SAT"):
				print("Valid solution.")
			elif line.startswith("UNSAT"):
				print("No valid solution. Returning")
				return
			else:
				solStr = line 

	ints = solStr.split(" ")[:N*dim]
	columnElems = []
	
	elem = K(0)
	for i in range(len(ints)):
		solComp = int(ints[i])
		exp = i % dim
		print i, exp, solComp
		#print solComp, i
		if solComp > 0:
			#print("Adding %s..."%str(w^exp))
			elem += w^exp
		if i%dim == dim - 1:
			print("Appending elem %s."%str(elem))
			columnElems.append(elem)
			elem = K(0)
	
	mat = sub.stack( matrix(k, columnElems) )
	mat = mat.augment(matrix(k, columnElems + [0]).transpose())
	return mat
	
def genMultiple(prefix, maxN):
	for i in range(maxN + 1):
		M = applySatSolution(prefix + "%d.txt"%i, A)
		lastCol = M[:, -1].list()
		generatePreSat(M,prefix + "%d.pre"%i, differ = lastCol)
	print("Saved all")
		
	
# Change last two columns of H, should be different:
A = H[:n-2, :n-2]
lastCol = vector( H[:,-2].list()[:-2] )
#print lastCol

generatePreSat(A, "sat%d.pre"%n, differ=lastCol)

		
	
