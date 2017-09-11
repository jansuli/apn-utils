from multiprocessing import cpu_count, Process, Manager, Event, Queue
from Queue import Empty
from numpy import array_split,array
import os

n = 8
K.<w> = GF(2^n, 'w', repr="log")
KSet = set(K.list())

dual = K.dual_basis()				# to get polynomial from QAM

resultDir = "QAMResults_n=%d"%n		# directory to save results in
if not os.path.isdir(resultDir):
	os.mkdir(resultDir)

nWorkers = max([4,cpu_count()])		# assume at least 4 cpu cores

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
		
def getPolynomFromQAM(qam):
	polStr = ""
	cf = Mtheta*qam*Mtheta.transpose()
	for j in range(n):
		for i in range(j):
			if cf[i,j] != 0:
				powerOfX = 2^i + 2^j
				polStr += "w^%s * x^%d + "%(str(cf[i,j]), powerOfX)
	return polStr[:-2]
			
# coefficient matrix for gold function	
Cf = matrix(K,n,n)
Cf [0,1] = Cf[1,0] = 1

H = M.transpose()*Cf*M

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
			complement = complement.difference(coset)
			
	return representatives

def getOriginalSetFunction(startMatrix):
	A = startMatrix
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
	Indices = [ind for ind in combinations[A.nrows()] if len(ind) > 0]	

	for ind in Indices:
		returnDict[tuple(ind)] = S(ind)

	def S(indices):
		return returnDict[indices]
		
	return S
			
#### Iterator as described in thesis ####

def getSetFn(oldFn, xi, pos):
	def newSetFn(indices):
		return oldFn(indices).intersection(set([xi + v for v in oldFn((pos,)+indices)]))
	return newSetFn	
					
def nextComponent(columnTilNow, setFn, nNeeded):
	newCol = copy(columnTilNow)
	pos = nNeeded + 1 - newCol.list().count(0) 
	domain = setFn((pos,))
	if len(domain) > 0:
		for xi in domain:
			newCol[pos] = xi			
			if pos < nNeeded:
				newSetFn = getSetFn(setFn, xi, pos)
				for col in nextComponent(newCol, newSetFn, nNeeded):
					yield col 
			else:
				yield newCol			
				
### Functions to be launched in parallel ###

def preSearchWorker(solutionQueue, startMatrix, adaptedSetFunction, enoughEvent, minSolutions):
	startColumn = matrix(K, startMatrix.nrows(), 1)
	solutionIterator = nextComponent(startColumn, adaptedSetFunction, startMatrix.nrows()-1)
	
	while True:
		try:
			solution = solutionIterator.next()
			solutionQueue.put(solution)
		
			if solutionQueue.qsize() >= minSolutions:
				enoughEvent.set()
		except StopIteration:
			print("preWorker %d done. Quitting."%os.getpid())
			break
				
def finalSearchWorkerInitiator(preSolutionQueue, startMatrix, nProcesses):
	while True:
		try:
			print("Initiating group of %d workers whilst having %d sub problems."%(nProcesses, preSolutionQueue.qsize()))
			preSol = preSolutionQueue.get(timeout=5)

			subMatrix = copy(startMatrix)
			N = startMatrix.nrows()

			subMatrix[:N-1, -1] = preSol
			subMatrix[-1, :N-1] = preSol.transpose()
			
			domainFunc = getOriginalSetFunction(subMatrix)
			adaptedDomainFunctions = partitionDomainFunc(domainFunc, nProcesses)
			
			workers = []
			for adaptedFx in adaptedDomainFunctions:		
				p = Process(target=finalSearchWorker, args=(subMatrix, adaptedFx))
				p.start()
				workers.append(p)
			for p in workers:
				p.join()
			
			print("All iterators stopped.")
		except Empty:
			print("No more sub problems.")
			break

def finalSearchWorker(subMatrix, adaptedFunc):
	print("%d starting work on \n%s\n"%(os.getpid(), subMatrix.str()))
	N = subMatrix.nrows()
	startColumn = matrix(K, N, 1)
	solutionIterator = nextComponent(startColumn, adaptedFunc, N-1)
	counter = 0
	while True:
		try:
			solution = solutionIterator.next()
			counter += 1
			
			B = matrix(K, N + 1, N + 1)
			B[:N,:N] = subMatrix
			B[-1, :N] = solution.transpose()
			B[:N, -1] = solution
			
			print("%d found QAM\n%s\n"%(os.getpid(), B.str()))
						
			with open(resultDir + "/res_%d.txt"%(os.getpid() + counter), "w") as f:
				f.write(B.str() +"\n" + getPolynomFromQAM(B))
		
		except StopIteration:
			print ("Process %d done, quitting."%os.getpid())
			break
		
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
	

if __name__ == "__main__":
	A2 = H[:n-2, :n-2]
	A1 = H[:n-1, :n-1]
	
	finWorkersPerSub = 4		# number of workers assigned to same subQAM	
	
	nPreWorkers = max( [1, floor(nWorkers/6)] )				# number of workers that change second last column
	nFinalWorkers = nWorkers - nPreWorkers
	nGroups = max(1,floor(nFinalWorkers/finWorkersPerSub))	# subQAMs to be checked in parallel
	
	S = getOriginalSetFunction(A2)
	adaptedSetFunctions = partitionDomainFunc(S, nPreWorkers)
	
	preSolutions = Queue()
	enoughEvent = Event()
	preWorkers = []
	print("Starting pre workers.")
	for adaptedFx in adaptedSetFunctions:
		p = Process( target=preSearchWorker, args=(preSolutions, A2, adaptedFx, enoughEvent, nFinalWorkers) )
		p.start()
		preWorkers.append(p)
	
	enoughEvent.wait()

	finalWorkers = []
	for i in range(nGroups):
		p = Process( target=finalSearchWorkerInitiator, args=(preSolutions, A1, finWorkersPerSub) )
		p.start()
		finalWorkers.append(p)
		
	for p in finalWorkers:
		p.join()
		
