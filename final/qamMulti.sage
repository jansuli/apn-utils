from multiprocessing import Process, Queue, cpu_count, Event, Value
from time import time
from os import getpid, path
from numpy import array_split

n = 6
K.<w> = GF(2^n, 'w',repr="log")
KSet = set(K.list())

dual = K.dual_basis()				# to get polynomial from QAM

resultDir = "QAMResults_%d"%n		# directory to save results in
if not path.isdir(resultDir):
	os.mkdir(resultDir)


nSplit = 3							# n_s in thesis
nCores = max( [nSplit, floor(cpu_count()/nSplit)*nSplit] ) 	# built for 24 cores
nGroups = nCores/nSplit				# n_p in thesis 

# Settings for adaptive Timeout
ADAPT_TIMEOUT = True		# wether or not timeout might change
INIT_STRATEGY = "MINIMUM" 	# "MINIMUM" or "MEAN", cf gaugeTimeout below
EASE_ON_FAIL  = True		# give workers more time to find QAMs when they fail
UPPER_FACTOR = 1.3			# if EASE_ON_FAIL: upper bound on timeout by multiplication with initial value
INCREMENTAL_FACTOR = 1.2	# if EASE_ON_FAIL: incremental factor on current timeout
INIT_FACTOR = 1.3			# factor to apply to value given by INIT_STRATEGY 
LOWER_BOUND = 100			# reset by gaugeTimeout
UPPER_BOUND = 1000			# reset by gaugeTimeout
WEIGHT = 1.1				# factor in calculation for mean of solution times

# determine combination indices beforehand to save time later on
combinations = {
	n-3 : Combinations(n-3).list(),
	n-2 : Combinations(n-2).list(),
	n-1 : Combinations(n-1).list(),
	n : Combinations(n).list(),	# nontrivial row combinations for qam checking
	}

# Transformation matrices M, Mtheta	
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
	'''Return polynomial (as string) to given n×n QAM.'''
	
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
	
print("Set up with initial QAM\n%s\nand polynomial %s."%(H.str(),getPolynomFromQAM(H)))
	
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
	
def getOriginalDomainFunction(submatrix):
	'''Takes a submatrix and returns the corresponding domain function S (as in thesis).'''
	
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
			
	# Possibly speed it up during search by saving results in dict
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
			S1parts = [set(arr) for arr in array_split(S1, m)]
						
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
	'''Apply solution given by iterator to submatrix and return now matrix.'''
	
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
			# MVR
			unassigned = sorted(unassigned, key = lambda v : len(domains[v]) )
			
			var = unassigned[0]
			for val in domains[var]:
				# Forwardcheck
				domainFn = getSetFn(domainFunc, val, var)
				newDomainDict = dict()
				for unVar in unassigned[1:]:
					domain = domainFn((unVar,))
					if domain == set():
						break
					newDomainDict[unVar] = domain
				else:
				# FC succeded. Continue ...
					assignment = copy(assigned)
					assignment[var] = val 
					for s in nextAssignement(assignment, newDomainDict, domainFn):
						if s != None: yield s
		else:
		# Last variable to assign
			assignment = copy(assigned)
			var = unassigned[0]
			for val in domains[var]:
				assignment[var] = val
				yield assignment

	solutionIterator = nextAssignement(assignment, domains, domainFunction)
	return solutionIterator

################################# Worker functions #####################

def preSearch(finalSolutionQueue, solutionIterator, startMatrix, initTimeout = None):
	''' Parent process to be launched with a Queue for final solutions, an iterator to begin with and a subQAM for fancy printing.'''
	
	pid = getpid()
	print("Process %d started."%pid)
	preSolutionQueue = Queue()
	
	preIterationStopped = False
	if initTimeout:		# has a timeout value been passed?
		timeoutObj = Value("d", initTimeout)	
	else: 
		TIMEOUT = None
		timeoutObj = None
	
	while True:
		if preSolutionQueue.qsize() > nSplit and preIterationStopped != True:
			print("\n--%s--\n-------------------------------------------------------"%pid)
			if timeoutObj:
				with timeoutObj.get_lock():
					TIMEOUT = timeoutObj.value
					 
			subprocesses = []
			subpids = []
			foundEvent = Event()											# Avoid killing productive workers
			for i in range(nSplit):
				finalDomainDict, finalSub = preSolutionQueue.get(timeout=5) # Get a job...
		
				p = Process( target=finalSearch, args=(foundEvent, finalDomainDict, finalSub, timeoutObj, finalSolutionQueue, pid) )
				p.start()
				subpids.append(str(p.pid))
				subprocesses.append(p)
			
			subpidStr = ",".join(subpids)
			print("Processes %s started final search."%subpidStr)
			subprocesses[0].join( timeout = TIMEOUT ) 						# join first process with timeout
			
			if foundEvent.is_set():
				# one subprocess found a solution, let all of them finish
				print("Waiting for all sub processes to finish")
				for p in subprocesses:
					p.join()
			else:
				downAdjust = False
				for p in subprocesses:
					if not p.is_alive():
						print("Sub process %d fininshed in time."%p.pid)
						downAdjust = True
					p.terminate()
				
				if ADAPT_TIMEOUT and timeoutObj and downAdjust:
					# one subprocess completed whole search before timeout -> timeout too low
					with timeoutObj.get_lock():
						newTIMEOUT = timeoutObj.value * 0.9
						timeoutObj.value = timeoutObj.value if newTIMEOUT < LOWER_BOUND else newTIMEOUT
						print("TIMEOUT for %d now is %d secs as it seemed to high."%(pid,timeoutObj.value))

				elif ADAPT_TIMEOUT and timeoutObj and EASE_ON_FAIL:		
					with timeoutObj.get_lock():
						newTIMEOUT =  floor(timeoutObj.value * INCREMENTAL_FACTOR)
						timeoutObj.value = timeoutObj.value if newTIMEOUT > UPPER_BOUND else newTIMEOUT
						print("TIMEOUT for %d now is %d secs as it seemed to harsh."%(pid,timeoutObj.value))
				if not downAdjust:
					print("All subprocesses timed out.")
		
		elif preIterationStopped != True:
			try:
				preSolution = solutionIterator.next()
				
				finalSubMatrix = extendSubQAM(A2, preSolution)
				finalS = getOriginalDomainFunction(finalSubMatrix)
				finalPartitionedFuncs = partitionDomainFunc(finalS, nSplit)
				
				for domainFunc in finalPartitionedFuncs:
					domainDict = dictFromSetFunction(domainFunc, finalSubMatrix.nrows())	# got to pass a dict instead of function
					preSolutionQueue.put( (domainDict, finalSubMatrix) )
					
			except StopIteration:
				preIterationStopped = True
		else:
			print("There are no more subproblems to be found or investigated. Process %d Stopping."%pid)
			break
			
def finalSearch(foundEv, domainDict, submatrix, timeoutObj, finalSolutionQueue, parentPID):
	'''Sub-process to find final QAM column.'''
	 
	pid = getpid()
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
				if ADAPT_TIMEOUT and timeoutObj:
					with timeoutObj.get_lock():
						oldTIMEOUT = timeoutObj.value
						newTIMEOUT = (q-1)/q * oldTIMEOUT + 1/q*floor(WEIGHT*(t_end - t_start))
						timeoutObj.value = newTIMEOUT if LOWER_BOUND <= newTIMEOUT <= UPPER_BOUND else oldTIMEOUT
						print("TIMEOUT for process %d now is %d secs."%(parentPID,timeoutObj.value))
					
			qam = extendSubQAM(submatrix, qamSolution)
			poly = getPolynomFromQAM(qam)
			finalSolutionQueue.put(poly)
			print("Subprocess %d of parent %d found a QAM:\n%s\n%s"%(pid,parentPID, qam.str(), poly))
			savePath = path.join(resultDir, "qam_%d_%d.txt"%(pid,solCounter))
			with open(savePath, "w") as f:
				f.write(qam.str() + "\n" + poly)
				print("Saved to '%s'."%savePath)
				
		except StopIteration:
			break
	
	print("Subprocess %d of parent %d stopping work."%(pid, parentPID))
	
def estimatorSearch(solIterator, firstSolTime):
	# launched by gaugeTimeout
	
	t_start = time()
	try:
		solIterator.next()
		t_end = time()
		print("Process %d found a solution in %f secs."%(getpid(), t_end - t_start))
		with firstSolTime.get_lock():
			firstSolTime.value = t_end - t_start
	except StopIteration:
		t_end = time()
		print("Process %d didn't find a solution within %f secs.."%(getpid(), t_end - t_start))
		with firstSolTime.get_lock():
			firstSolTime.value = t_start - t_end	# negative sign for distinction
		
		
def gaugeTimeout():
	'''Estimate timeout by investigating gold QAM.'''
	
	print("Trying to append to gold QAM to gauge timeout.")
	global LOWER_BOUND, UPPER_BOUND
	
	A = H[:n-1, :n-1]
	S = getOriginalDomainFunction(A)
	domainFuncs = partitionDomainFunc(S, nSplit)
	estimators = []
	timesV = []
	for fx in domainFuncs:
		sol = setUpIterator(n-1, fx)
		solTime = Value('d')
		p = Process( target=estimatorSearch, args=(sol,solTime) )
		p.start()
		estimators.append(p)
		timesV.append(solTime)
	for p in estimators:
		p.join()
	
	print("Estimations done.")
	times = []
	for sT in timesV:
		with sT.get_lock():
			times.append(sT.value)
			
	solTimes = [t for t in times if t >= 0 ]
	unSolTimes = [-t for t in times if t < 0]
	
	if len(solTimes) > 0:
		if INIT_STRATEGY == "MINIMUM":
			timeout = max([1,floor(min(solTimes)*INIT_FACTOR)])
	
		elif INIT_STRATEGY == "MEAN":
			timeout = max([1,floor(mean(solTimes)*INIT_FACTOR)])
	else:
		timeout = floor(min(unSolTimes))
		
	LOWER_BOUND = timeout
	UPPER_BOUND = ceil(timeout*UPPER_FACTOR)
	
	print("Initial Timeout set to %d secs by %s strategy."%(timeout, INIT_STRATEGY))
	
	return timeout
	
def saveSols(solQueue):
	'''Save all solutions currently in the given queue as comma-seperated list in .txt file.'''
	 
	print("Saving solutions...")
	sols = ""
	while solQueue.qsize() > 0:
		poly = solQueue.get()
		sols = sols + poly + ",\n"
	if sols != "":
		with open(path.join(resultDir, "resList.txt"), "w") as f:
			f.write(sols)
	print("Done.")
		
if __name__ == "__main__":
	try:
		initialTimeout = gaugeTimeout()
		print(LOWER_BOUND, UPPER_BOUND)
		# Set up pre-Search preliminaries
		A2 = H[:n-2, :n-2]
		S2 = getOriginalDomainFunction(A2)
		partitionedFuncs = partitionDomainFunc(S2, nGroups)
		organizers = []
		qamSolutions = Queue()
		
		for j in range(nGroups):
			print("Setting up group %d of workers."%j)
			# each process gets its own iterator
			preIterator = setUpIterator(A2.nrows(), partitionedFuncs[j] )
			p = Process( target=preSearch, args=(qamSolutions, preIterator, A2, initialTimeout) )
			p.start()
			organizers.append(p)
		
		for p in organizers:
			p.join()
		
		saveSols(qamSolutions)
	except KeyboardInterrupt:
		saveSols(qamSolutions)
		
			
