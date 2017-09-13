from multiprocessing import Process, Queue, cpu_count, Event, Value
from time import time
from os import getpid, path
from numpy import array_split

n = 8
K.<w> = GF(2^n, 'w',repr="log")
KSet = set(K.list())

dual = K.dual_basis()				# to get polynomial from QAM

resultDir = "QAMResults_%d"%n		# directory to save results in
if not path.isdir(resultDir):
	os.mkdir(resultDir)


nSplit = 3
nCores = max( [nSplit, floor(cpu_count()/nSplit)*nSplit] ) 	# built for 24 cores
nGroups = nCores/nSplit

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
		
Mtheta = matrix(K,n,n)
for i in range(n):
	for j in range(n):
		Mtheta[i,j] = dual[j]^(2^i)
		
# coefficient matrix for gold function	
Cf = matrix(K,n,n)
Cf[0,1] = Cf[1,0] = 1

H = M.transpose()*Cf*M

def getPolynomFromQAM(qam):
	polStr = ""
	cf = Mtheta*qam*Mtheta.transpose()
	for j in range(n-1,-1, -1):
		for i in range(j-1,-1, -1):
			if cf[i,j] != 0:
				powerOfX = 2^i + 2^j
				if cf[i,j] != 2^n-1:
					polStr += "g^%s*x^%d + "%(str(cf[i,j]), powerOfX)
				else:
					polStr += "x^%d + "%(powerOfX)
	return polStr[:-3]
	
print(getPolynomFromQAM(H))
	
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

def preSearch(finalSolutionQueue, solutionIterator, startMatrix, initTimeout = None):
	pid = getpid()
	print("Process %d started."%pid)
	preSolutionQueue = Queue()
	
	preIterationStopped = False
	if initTimeout:
		timeoutObj = Value("d", initTimeout)									# has a timeout value been passed?
	else: 
		TIMEOUT = None
		timeoutObj = None
	
	timeoutCounter = 0
	while True:
		if preSolutionQueue.qsize() > 1 and preIterationStopped != True:
			if timeoutObj:
				with timeoutObj.get_lock():
					TIMEOUT = timeoutObj.value
					 
			finalDomainDict, finalSub = preSolutionQueue.get(timeout=5) 		# Get a job...
		
			foundEvent = Event()												# Avoid killing productive workers
			p = Process( target=finalSearch, args=(foundEvent, finalDomainDict, finalSub, timeoutObj, finalSolutionQueue, pid) )
			p.start()
			t_start = time()
			p.join( timeout=TIMEOUT )
			if p.is_alive():
				print("Subprocess %d of Parent %d timed out."%(p.pid, pid))
				if foundEvent.is_set():
					print("But it found a solution. So we will let it run.")
					timeoutCounter = 0
					p.join()
				else:
					p.terminate()
					timeoutCounter += 1
					if timeoutCounter == nSplit:
						with timeoutObj.get_lock():
							timeoutObj.value = timeoutObj.value * 1.25
							print("TIMEOUT for %d now is %d secs as it seemed to harsh."%(pid,timeoutObj.value))
						timeoutCounter = 0
			else:
				print("Subprocess %d of Parent %d stopped within timeout."%(p.pid, pid))
				t_end = time()
				tDelta = 0.9* (t_end - t_start)
				with timeoutObj.get_lock():
					timeoutObj.value = tDelta
					print("TIMEOUT for %d lowered from %d to %d."%(pid,TIMEOUT, tDelta))
						
		
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
			
def finalSearch(foundEv, domainDict, submatrix, timeoutObj, finalSolutionQueue, parentPID):
	pid = getpid()
	print("%d started looking for final column."%pid)
	def domainFunc(ind):
		return domainDict[ind]
	
	solutionIterator = setUpIterator(submatrix.nrows(), domainFunc)
	solCounter = 0
	t_start = time()
	solCounter = 0
	while True:
		try:
			qamSolution = solutionIterator.next()
			solCounter += 1
			foundEv.set()
			
			if solCounter == 1:
				t_end = time()
				q = finalSolutionQueue.qsize() + 2
				with timeoutObj.get_lock():
					oldTIMEOUT = timeoutObj.value 
					newTIMEOUT = (q-1)/q * oldTIMEOUT + 1/q*floor(1.1*(t_end - t_start))
					timeoutObj.value = newTIMEOUT
				print("TIMEOUT for %d now is %d secs."%(parentPID,newTIMEOUT))
					
			qam = extendSubQAM(submatrix, qamSolution)
			poly = getPolynomFromQAM(qam)
			finalSolutionQueue.put(qam)
			print("%d found a QAM:\n%s\n%s"%(pid, qam.str(), poly))
			savePath = path.join(resultDir, "qam_%d_%d.txt"%(pid,solCounter))
			with open(savePath, "w") as f:
				f.write(qam.str() + "\n" + poly)
				print("Saved to '%s'."%savePath)
				
		except StopIteration:
			break
	
	print("%d stopping work."%pid)
	
def estimatorSearch(solIterator, firstSolTime):
	t_start = time()
	try:
		solIterator.next()
		t_end = time()
		print("Process %d found a solution in %f secs."%(getpid(), t_end - t_start))
		with firstSolTime.get_lock():
			firstSolTime.value += t_end - t_start
	except StopIteration:
		t_end = time()
		print("Process %d didn't find a solution within %f secs.."%(getpid(), t_end - t_start))
		with firstSolTime.get_lock():
			firstSolTime.value += t_start - t_end	# negative sign for distinction
		
		
def gaugeTimeout():
	print("Trying to append to gold QAM to gauge timeout.")
	
	A = H[:n-1, :n-1]
	S = getOriginalDomainFunction(A)
	domainFuncs = partitionDomainFunc(S, nSplit)
	estimators = []
	times = []
	for fx in domainFuncs:
		sol = setUpIterator(n-1, fx)
		solTime = Value('d')
		p = Process( target=estimatorSearch, args=(sol,solTime) )
		p.start()
		estimators.append(p)
		times.append(solTime)
	for p in estimators:
		p.join()
		
	print("Estimations done.")
	solMean = 0
	divisorSol = 0
	unSolMean = 0
	divisorUnSol = 0
	for t in times:
		if t.value > 0:
			solMean += t.value
			divisorSol += 1
		else:
			unSolMean += -t.value
			divisorUnSol += 1
			
	if divisorSol != 0 and divisorUnSol != 0: 
		solMean = solMean/divisorSol
		unSolMean = unSolMean/divisorUnSol
		
	if divisorSol != 0 and divisorUnSol != 0:
		factor = (unSolMean/solMean - 1)/2 + 1
		timeout = solMean * factor
	elif divisorSol != 0 and divisorUnSol == 0:
		timeout = 3*solMean
	elif divisorSol == 0 and divisorUnSol != 0:
		timeout = unSolMean
	timeout = floor(timeout)				# milliseconds
	print("Mean time for finding first solution: %f secs."%solMean)
	print("Mean time for doing complete search: %f secs."%unSolMean)
	print("Initial Timeout set to %d secs."%timeout)
	
	return timeout
		
if __name__ == "__main__":
	initialTimeout = gaugeTimeout()
	
	# Set up pre-Search preliminaries
	A2 = H[:n-2, :n-2]
	S2 = getOriginalDomainFunction(A2)
	partitionedFuncs = partitionDomainFunc(S2, nCores)
	organizers = []
	qamSolutions = Queue()
	for j in range(nGroups):
		print("Setting up group %d of workers."%j)
		for i in range(nSplit):
			# each process gets its own iterator
			preIterator = setUpIterator(A2.nrows(), partitionedFuncs[i+j] )
			p = Process( target=preSearch, args=(qamSolutions, preIterator, A2, initialTimeout) )
			p.start()
			organizers.append(p)
	
