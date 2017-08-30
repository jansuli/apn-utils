from multiprocessing import cpu_count, Process, Manager
from numpy import array_split,array

n = 9
K.<w> = GF(2^n, 'w',repr="log")
KSet = set(K.list())

nWorkers = cpu_count()

# determine combination indices beforehand to save time later on
combinations = {
	n-3 : Combinations(n-3).list(),
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

def rowSpan(rowVector):
	'''Takes a row vector and outputs the subspace spaned by its elements as a set.'''
	
	N = len(rowVector)
	print N
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
			complement = complement.difference(coset)
			
	return representatives

def getOriginalSetFunction(startMatrix):
	A = startMatrix
	# search space reduction
	print("First position.")
	S0 = cosetRepr(A.rows()[0])
	print S0
	S1 = set()
	print("Second position.")
	Vred = rowSpan(A.rows()[1][1:])
	for v in cosetRepr(A.rows()[1]):
		print("Repr %s."%str(v))
		S1 = S1.union(set( v + a for a in Vred))

	print("Defininig")
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
			
	print("Saving in dict.")
	# Possibly speed it up during search
	returnDict = dict()
	Indices = [ind for ind in combinations[A.nrows()] if len(ind) > 0]	

	for ind in Indices:
		returnDict[tuple(ind)] = S(ind)

	print("Done.")
	def S(indices):
		return returnDict[indices]
		
	return S
			
####

def getSetFn(oldFn, xi, pos):
	def newSetFn(indices):
		indices = tuple(sorted(indices))
		return oldFn(indices).intersection(set([xi + v for v in oldFn(tuple([pos])+indices)]))
	return newSetFn
	
def nextAssignement(variables, assigned, domains, domainFunc, nRows):
	unassigned = [v for v in variables if v not in assigned]

	if len(assigned) < nRows-1:
	# We still need to assign at least one variable
		
		# Sorting by domain size for M(ost) R(estricted) V(ariable) heuristic
		unassigned = sorted(unassigned, key = lambda v : len(domains[v]) )
		
		# Try to assign values to unassigned variables
		for var in unassigned:
			unassignedLeft = [unVar for unVar in unassigned if unVar != var]
			for val in domains[var]:
				
				# Forwardcheck : no unassigned variable should have an empty domain
				domainFn = getSetFn(domainFunc, val, var)
				newDomainDict = dict()
				for unVar in unassignedLeft:
					domain = domainFn((unVar,))
					if domain == set():
						# domain empty
						break
					newDomainDict[unVar] = domain
				else:
					# Forwardcheck succeeded! Recursion
					assignment = copy(assigned)
					assignment[var] = val 
					for s in  nextAssignement(variables, assignment, newDomainDict, domainFn, nRows):
						yield s
	else:
		# Only one unassigned variable with non-empty domain (ensured beforhand by FC)
		assignment = copy(assigned)
		var = unassigned[0]
		for val in domains[var]:
			assignment[var] = val
			yield assignment

def searchWorker(solutionList, startMatrix, adaptedDomainFunc):
	assignment = dict()
	variables = set(range(startMatrix.nrows()))
	domains = dict()
	for v in variables:
		domains[v] = adaptedDomainFunc( (v,) )
		
	print variables,domains

	sol = nextAssignement(variables, assignment, domains, adaptedDomainFunc, startMatrix.nrows())
	
	while True: print sol.next()
	#while len(solutionList) < nWorkers:
		#solutionList.append(sol.next())
		
### Partition search space
A2 = H[:n-2, :n-2]
S = getOriginalSetFunction(A2)

S0 = list(S((0,)))
S1 = list(S((1,)))

adaptedSetFunctions = []
if len(S0)/nWorkers >= 1:
	# more (or equally as many) options than workers
	S0parts = [set(arr) for arr in array_split(S0, nWorkers)]
	print("For %d workers splitting only S0 in %d parts."%(nWorkers, nWorkers))	
	
	def getNewSetFunc(part):
		def newFunc(indices):
			if indices == (0,):
				return part
			else:
				return S(indices)
		return newFunc
		
	for part in S0parts:
		adaptedSetFunctions.append(getNewSetFunc(part))
		
else:
	def getNewSetFunc(part0, part1):
			def newFunc(indices):
				if indices == (0,):
					return part0
				elif indices == (1,):
					return part1
				else:
					return S(indices)
			return newFunc
			
	r = nWorkers % len(S0)
	if r == 0:
		m = nWorkers / len(S0)
		S0parts = [set([elem]) for elem in S0]
		S1parts = [set(arr) for arr in array_split(S0, m)]
		print("For %d splitting S0 in %d and S1 in %d parts."%(nWorkers, nWorkers, secondSplit))
		
		for part0 in S0parts:
			for part1 in S1parts:
				adaptedSetFunctions.append(getNewSetFunc(part0, part1))
	
	else:
		m = int(floor(nWorkers/len(S0)))
		print("For %d Workers we get:"%nWorkers)
		print("%d processes for which S1 is split into %d parts …"%((len(S0)-r), m))
		S0 = list(S0)
		print len(S0[r:])
		for elem in S0[r:]:
			for arr in array_split(S1,m):
				adaptedSetFunctions.append(getNewSetFunc(set([elem]), set(arr)))	
		print("… and %d processes for which S1 is split into %d parts."%(r, (m+1)))
		print len(S0[:r])
		for elem in S0[:r]:
			for arr in array_split(S1, m+1):
				adaptedSetFunctions.append(getNewSetFunc(set([elem]), set(arr)))
		
m = Manager()
solutions = m.list()
workers = []

if __name__ == "__main__":

	for partitionFunc in adaptedSetFunctions:
		p = Process(target = searchWorker, args = (solutions, A2, partitionFunc))
		p.start()
		workers.append(p)
		
