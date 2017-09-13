from multiprocessing import Process, Queue, cpu_count, Event, Value
from time import time
from os import getpid
from numpy import array_split

n = 7
K.<w> = GF(2^n, 'w',repr="log")
KSet = set(K.list())

nSplit = 3
nCores = max( [nSplit, floor(cpu_count()/nSplit)*nSplit] ) # built for 24 cores
nGroups = nCores/nSplit

TIMEOUT = 50000 # in milliseconds

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
Cf[0,1] = Cf[1,0] = 1

H = M.transpose()*Cf*M

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

def rowSpan(rowVector):
	'''Takes a row vector and outputs the subspace spaned by its elements as a set.'''
	
	N = len(rowVector)
	span = set()
	for ind in combinations[N]:
		newElem = sum([rowVector[i] for i in ind])
		span.add( newElem )
	return span
	
def cosetRepr(rowVector):
	'''Takes a row vector and returns a set of representative elements of nontrivial cosets of its span.'''
	 
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
	
def getOriginalDomainFunction(submatrix, returnDict = False):
	N = submatrix.nrows()

	# search space reduction
	S0 = cosetRepr(submatrix.rows()[0])
	S1 = set()
	Vred = rowSpan(submatrix.rows()[1][1:])
	for v in cosetRepr(submatrix.rows()[1]):
		S1 = S1.union(set( v + a for a in Vred))

	def S(indices):
		if len(indices) == 1 :
			index = indices[0]
			if index == 0:
				return S0
			elif index == 1:
				return S1
			else:
				return KSet.difference( rowSpan(submatrix.rows()[index]) )
		else:
			return KSet.difference( rowSpan( sum(submatrix[indices, :]) ) )
			
	# Possibly speed it up during search
	returnDict = dict()
	Indices = [ind for ind in combinations[N] if len(ind) > 0]	

	for ind in Indices:
		returnDict[tuple(ind)] = S(ind)

	def S(indices):
		return returnDict[indices]

	return S
	
## avoid pickling functions
def dictFromSetFunction(setFunction, indexKey):
	returnDict = dict()
	for ind in combinations[indexKey]:
		if len(ind) > 0:
			ind = tuple(ind)
			returnDict[ind] = setFunction(ind)
	return returnDict
	
# Split a given domain function
def partitionDomainFunc(originalSetFunc, nProcesses):
	Sf = originalSetFunc
	S0 = list(Sf((0,)))
	S1 = list(Sf((1,)))

	adaptedSetFuncs = []
	if len(S0)/nProcesses >= 1:
		# more (or equally as many) options than workers
		S0parts = [set(arr) for arr in array_split(S0, nProcesses)]	
		
		def getNewSetFunc(part):
			def newFunc(indices):
				if indices == (0,):
					return part
				else:
					return Sf(indices)
			return newFunc
			
		for part in S0parts:
			adaptedSetFuncs.append(getNewSetFunc(part))
			
	else:
		def getNewSetFunc(part0, part1):
				def newFunc(indices):
					if indices == (0,):
						return part0
					elif indices == (1,):
						return part1
					else:
						return Sf(indices)
				return newFunc
				
		r = nProcesses % len(S0)
		if r == 0:
			m = nProcesses / len(S0)
			S0parts = [set([elem]) for elem in S0]
			S1parts = [set(arr) for arr in array_split(S0, m)]
						
			for part0 in S0parts:
				for part1 in S1parts:
					adaptedSetFuncs.append(getNewSetFunc(part0, part1))
		else:
			m = int(floor(nProcesses/len(S0)))
			S0 = list(S0)
			for elem in S0[r:]:
				for arr in array_split(S1,m):
					adaptedSetFuncs.append(getNewSetFunc(set([elem]), set(arr)))	
			for elem in S0[:r]:
				for arr in array_split(S1, m+1):
					adaptedSetFuncs.append(getNewSetFunc(set([elem]), set(arr)))
					
	return adaptedSetFuncs
	
def extendSubQAM(submatrix, solution):
	N = submatrix.nrows()
	B = matrix(K, N+1)
	B[:N, :N] = submatrix
	for pos, value in solution.iteritems():
		B[-1, pos] = B[pos, -1] = value
	return B
	
def setUpIterator(nRows, domainFunction):
	#### MVR+FC Iterator 

	assignment = dict()
	variables = set(range(nRows))
	domains = dict()

	def getSetFn(oldFn, xi, pos):
		def newSetFn(indices):
			indices = tuple(sorted(indices))
			return oldFn(indices).intersection(set([xi + v for v in oldFn((pos,)+indices)]))
		return newSetFn
			
	for v in variables:
			domains[v] = domainFunction( (v,) )

	def nextAssignement(assigned, domains, domainFunc):
		unassigned = [v for v in variables if v not in assigned]

		if len(unassigned) > 1 :
			# mvr:
			unassigned = sorted(unassigned, key = lambda v : len(domains[v]) )
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
			var = unassigned[0]
			for val in domains[var]:
				assignment[var] = val
				yield assignment

	solutionIterator = nextAssignement(assignment, domains, domainFunction)
	return solutionIterator

################################# Worker functions #####################

def preSearch(preSolutionQueue, solutionIterator, startMatrix, timeoutObj = None):
	pid = getpid()
	print("Process %d started."%pid)
	preIterationStopped = False
	 
	while True:
		if preSolutionQueue.qsize() > 1 and preIterationStopped != True:
			finalDomainDict, finalSub = preSolutionQueue.get(timeout=5) 
		
			foundEvent = Event()
			p = Process( target=finalSearch, args=(foundEvent, finalDomainDict, finalSub) )
			p.start()
			p.join( timeout=timeoutObj )
			if p.is_alive():
				print("Subprocess timed out.")
				if foundEvent.is_set():
					print("But it found a solution. So we will let it run.")
					p.join()
				else:
					print("And no solution. Killing.")
					p.terminate()
			else:
				print("Subprocess stopped within timeout.")
		
		elif preIterationStopped != True:
			try:
				preSolution = solutionIterator.next()
				
				finalSubMatrix = extendSubQAM(A2, preSolution)
				finalS = getOriginalDomainFunction(finalSubMatrix)
				finalPartitionedFuncs = partitionDomainFunc(finalS, nSplit)
				
				for domainFunc in finalPartitionedFuncs:
					domainDict = dictFromSetFunction(domainFunc, finalSubMatrix.nrows())
					preSolutionQueue.put( (domainDict, finalSubMatrix) )
					
			except StopIteration:
				preIterationStopped = True
		else:
			print("There are no more subproblems to be found or investigated. %d Stopping."%pid)
			break
			
def finalSearch(foundEv, domainDict, submatrix):
	pid = getpid()
	print("%d started looking for final column."%pid)
	def domainFunc(ind):
		return domainDict[ind]
	
	solutionIterator = setUpIterator(submatrix.nrows(), domainFunc)
	while True:
		try:
			qamSolution = solutionIterator.next()
			foundEv.set()
			qam = extendSubQAM(submatrix, qamSolution)
			print("%d found a QAM:\n%s\n."%(pid, qam.str()))
			checkQAM(qam)
		except StopIteration:
			break
	
	print("%d stopping work."%pid)
	
def gaugeTimeout(timeoutObject):
	print("Trying to append to gold QAM to gauge timeout.")
	solutionTime = 0
	divisor = 0
	A = H[:n-1, :n-1]
	S = getOriginalDomainFunction(A)
	domainFuncs = partitionDomainFunc(S, nSplit)
	for fx in domainFuncs:
		sol = setUpIterator(n-1, fx)
		t_start = time()
		try:
			sol.next()
			t_end = time()
			solutionTime += (t_end - t_start)
			divisor += 1
		except StopIteration:
			t_end = time()
			solutionTime += (t_end - t_start)
			divisor += 1
			
	if divisor != 0:
		timeout = floor(solutionTime/divisor)*2000
	else:
		timeout = 0
		
	with timeoutObject.get_lock():
		timeoutObject = timeout
		print(timeout)
		
if __name__ == "__main__":
	timeoutObject = Value('d')
	#gaugeTimeout(timeoutObject)
	# Set up pre-Search preliminaries
	A2 = H[:n-2, :n-2]
	S2 = getOriginalDomainFunction(A2)
	partitionedFuncs = partitionDomainFunc(S2, nCores)
	organizers = []
	for j in range(nGroups):
		preSolutions = Queue()
		for i in range(nSplit):
			print i
			# each process gets its own iterator
			preIterator = setUpIterator(A2.nrows(), partitionedFuncs[i+j] )
			p = Process( target=preSearch, args=(preSolutions, preIterator, A2) )
			p.start()
			organizers.append(p)
	
