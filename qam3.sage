import pickle

n = 7
K.<w> = GF(2^n, 'w' )#,repr="log")
KSet = set(K.list())
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
	N = len(row)
	for ind in combinations[N]:
		ret.add( K(sum([row[i] for i in ind])) )
	return ret	
	
def f(x):
	return x^3
	
Cf = matrix(K,n,n)
Cf [0,1] = Cf[1,0] = 1

H = M.transpose()*Cf*M

def rowSpan(rowVector):
	'''Takes a row vector and outputs the subspace spaned by its elements as a set.'''
	
	N = len(rowVector)
	span = set()
	for ind in combinations[N]:
		newElem = sum([rowVector[i] for i in ind])
		span.add( K(newElem) )
	return span
		
def cosetRepr(rowVector):
	'''Takes a row vector and returns a set of representative elements of cosets translating the row space.'''
	 
	V = rowSpan(rowVector)
	complement = KSet.difference(V)
	representatives = set()
	build = set()
	
	V = list(V)
	a = choice(V)
	while a == 0:
		a = choice(V)
	for w in complement:
		if w + a not in build:
			coset = set( [v+w for v in V] )
			build = build.union(coset)
			representatives.add(w)
	return representatives

		
def differConstraint(vec):
	'''Takes an n-1 vector of F_{2^n} elements and returns a boolean function as CNF expression in (n-1)*n variables ensuring that the binary vectors differ.'''
	expression = "And( "
	count = 0
	for elem in vec:
		subexpression = "Or( "
		binary = vector(elem)
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
	
	# Add domain constraints using S
	# search space reduction
	S0 = cosetRepr(A.rows()[0])
	S1 = set()
	Vred = rowSpan(A.rows()[1][1:])
	for v in cosetRepr(A.rows()[1]):
		S1 = S1.union(set( v + a for a in Vred))

	def S(indices):
		if len(indices) == 1 :
			index = indices[0]
			if index == 0:
				return S0
			elif index == 1:
				return S1
			else:
				return KSet.difference( rowSpan(A.rows()[index]) )
		else:
			return KSet.difference( rowSpan( sum(A[indices, :]) ) )
			
	for i in range(N):
		offset = i*dim + 1
		
		formatting = ["x%d"%j for j in range(offset, offset+dim)]
		expr = notInSetConstraint(KSet.difference(S((i,)))) 
		constraintExpression = expr.format(*formatting)
		expressionList.append(constraintExpression)
		
	expression1 = "And("+ ",".join(expressionList)+")"  # everything up till now is already in CNF and doesn't need to be transformed
	expressionList = []
	
	sumIndices = [ind for ind in combinations[N] if len(ind) > 1]
	for ind in sumIndices:
		formatting = ["Xor( " for i in range(dim)] 	# put n xor's into constraintFunctions  
		
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
		
	
# Change last two columns of H, should be different:
A = H[:n-1, :n-1]
lastCol = vector( H[:,-1].list()[:-1] )

generatePreSat(A, "_sat%d.pre"%n) # differ=lastCol)
