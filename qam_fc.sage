n = 7
K.<w> = GF(2^n, 'w',repr="log")
KSet = set(K.list())

# determine combination indices beforehand to save time later on
combinations = {
	n-2 : Combinations(n-2).list(),
	n-1 : Combinations(n-1).list(),
	n : Combinations(n).list(),	# nontrivial row combinations for qam checking
	}
	
M = matrix(K, n, n)

for i in range(n):
	for j in range(n):
		M[i,j] = w^(j*2^i)
		
# coefficient matrix for gold function	
Cf = matrix(K,n,n)
Cf [0,1] = Cf[1,0] = 1

H = M.transpose()*Cf*M

A = H[:n-1, :n-1]

def rowSpan(rowVector):
	'''Takes a row vector and outputs the subspace spaned by its elements as a set.'''
	
	N = len(rowVector)
	span = set()
	for ind in combinations[N]:
		newElem = sum([rowVector[i] for i in ind])
		span.add( newElem )
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
		
# Possibly speed it up during search
returnDict = dict()
Indices = [ind for ind in combinations[n-1] if len(ind) > 0]	

for ind in Indices:
	returnDict[tuple(ind)] = S(ind)

def S(indices):
	return returnDict[indices]
		
#### MVR+FC Iterator 

assignment = dict()
variables = set(range(A.nrows()))
domains = dict()

def getSetFn(oldFn, xi, pos):
	def newSetFn(indices):
		indices = tuple(sorted(indices))
		return oldFn(indices).intersection(set([xi + v for v in oldFn((pos,)+indices)]))
	return newSetFn
		
for v in variables:
		domains[v] = S( (v,) )

def nextAssignement(assigned, domains, domainFunc):
	unassigned = [v for v in variables if v not in assigned]

	if len(assigned) < n-2:
		# mvr:
		#unassigned = sorted(unassigned, key = lambda v : len(domains[v]) )
		var = unassigned[0]
		unassignedLeft = [unVar for unVar in unassigned if unVar != var]
		for val in domains[var]:
			# forwardcheck
			domainFn = getSetFn(domainFunc, val, var)
			newDomainDict = dict()
			for unVar in unassignedLeft:
				domain = domainFn((unVar,))
				if domain == set():
					#print("Forward Check failed. Backtracking.")
					break
				newDomainDict[unVar] = domain
			else:
				assignment = copy(assigned)
				assignment[var] = val 
				for s in nextAssignement(assignment, newDomainDict, domainFn):
					if s != None: yield s
	else:
		assignment = copy(assigned)
		for val in domains[unassigned[0]]:
			assignment[unassigned[0]] = val
			yield assignment

sol = nextAssignement(assignment, domains, S)

### Original Iterator
			
def getSetFnOr(oldFn, xi, pos):
	def newSetFn(indices):
		return oldFn(indices).intersection(set([xi + v for v in oldFn((pos,)+indices)]))
	return newSetFn	
					
def nextComponent(columnTilNow, setFn):
	newCol = copy(columnTilNow)
	pos = n - (newCol.list().count(0) + 1)
	domain = setFn(tuple([pos]))
	if len(domain) > 0:
		for xi in domain:
			newCol[pos] = xi			
			if pos < n-2:
				newSetFn = getSetFnOr(setFn, xi, pos)
				for col in nextComponent(newCol, newSetFn):
					if col != None: yield col 
			else:
				yield newCol
	else:
		yield

startColumn = matrix(K, n-1, 1)		
sol2 = nextComponent(startColumn, S)
		

def applySol(assignment):
	B = copy(H)
	for pos, val in assignment.iteritems():
		B[pos, -1] = B[-1, pos] = val
	return B
		
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
		
def checkRowRank(row):
	testMatrix = matrix( GF(2), n,0)
	for i in range(row.ncols()):
		testMatrix = testMatrix.augment(matrix(vector(row[0,i])).transpose())
	#print testMatrix.str()
	return testMatrix.rank()
		

