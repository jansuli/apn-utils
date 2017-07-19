#from sympy import var as syVar
from sympy.logic.boolalg import to_cnf, conjuncts
from tqdm import tqdm
import pickle

n = 5
K.<w> = GF(2^n, 'w' )#,repr="log")
k = K.list()
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
	
def cosetRepr(space, subspace):
	subspace = list(subspace)
	num = 2^(log(len(space), base=2) - log(len(subspace) , base=2))
	reps = [K(0)]
	cosetsUnion = set(subspace)
	i = 0
	v = subspace[1]
	while len(reps) < num:
		a = space[i]
		if a != 0:
			if not a+v in cosetsUnion:
				reps.append(a)
				newCoset =  set([a + w for w in subspace]) 
				cosetsUnion = cosetsUnion.union(newCoset)
				
		i += 1
	return reps
	
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
	
def notInSetConstraint(outSet, index):
	expression = ""
	vectorIndices = list(range(index*n + 1, index*n+n + 1))
	for elem in outSet:
		subexpression = ""
		vec = vector(elem)
		count = 0
		for entry in vec:
			if entry == 1:
				subexpression += "-{%d} "%count
			else:
				subexpression += "{%d} "%count
			count += 1
		expression += (subexpression + " 0\n").format(*vectorIndices)
	#expression = expression[:-1]	#remove last \n
	print ("Added %d clauses."%len(outSet))
	return expression, len(outSet)
	
def notInCosetProduct(reprTuple, subspace, vectorIndices):
	cosets = []
	for rep in reprTuple:
		cosets.append( Set( [rep + w for w in subspace] ) )
	fieldTuples = cosets[0].cartesian_product(*cosets[1:])
	
	expression = ""
	clauses = 0
	for fieldTuple in fieldTuples:
		longVec = []
		subexpression = ""
		count = 0
		for elem in fieldTuple:
			longVec += vector(elem).list()
		#print longVec, len(longVec)
		for entry in longVec:
			if entry == 1:
				subexpression += "-{%d} "%count
			else:
				subexpression += "{%d} "%count
			count += 1
		expression += (subexpression + " 0\n").format(*vectorIndices)
		clauses += 1
	#expression = expression[:-1]
	
	return expression, clauses
	
def sumNotInSet(outSet, ind):
	vectorIndices = []
	for i in ind:
		vectorIndices += list(range(i*n+1, i*n+n+1))
	#print vectorIndices
	
	reprTuple = cosetRepr(k, outSet)
	
	lenTuple = len(ind)
	variations = Tuples(reprTuple, lenTuple)
	expression = ""
	clauses = 0
	for var in variations:
		if lenTuple % 2 == 0:
			for rep in var:
				if var.count(rep) % 2 == 1:
					break
			else:
				# sum of representatives is 0 => sum in outSet => f would evaluate to 0
				print("not in coset product constraint")
				print("not in coset product constraint")
				print (ind)
				print (var)
				print ("\n")
				
				subexpression, subclauses = notInCosetProduct(var, outSet, vectorIndices)
				expression += subexpression
				clauses += subclauses
		elif lenTuple % 2 == 1 and 0 in var:
			cp = var[:]
			cp.remove(0)
			for rep in cp:
				if cp.count(rep) % 2 == 1:
					break
			else:
				# sum of representatives is 0 => sum in outSet => f would evaluate to 0
				print(var)
				print("not in coset product constraint")
				print (ind)
				print (var)
				print ("\n")
				
				subexpression, subclauses = notInCosetProduct(var, outSet, vectorIndices)
				expression += subexpression
				clauses += subclauses
	
	return expression, clauses
			 
def generateSatString(submatrix):
	N = submatrix.nrows()
	d = submatrix.base_ring().degree()
	expression = ""
	# domain constraints:
	rows = submatrix.rows()
	clauses = 0
	for i in range(0, N):
		rowSpace = rowSpan(rows[i])
		print len(rowSpace)
		subexpression, subclauses = notInSetConstraint(rowSpace, i)
		expression += subexpression
		clauses += subclauses
	if True:	
		sumIndices = [ind for ind in combinations[N] if len(ind) > 1]
		for ind in sumIndices:
			rowSpace = rowSpan(sum([rows[i] for i in ind]))
			subexpression, subclauses = sumNotInSet(rowSpace, ind)
			expression += subexpression
			clauses += subclauses
		
	return "p cnf %d %d \n"%(d*N, clauses) + expression[:-1]

A = H[:n-1, :n-1]
s = generateSatString(A)

with open("directsat%d.cnf"%n, "w") as f:
	f.write(s)
	
	
	
	
	


