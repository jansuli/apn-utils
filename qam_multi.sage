from multiprocessing import cpu_count, Process, Manager, Event
from numpy import array_split,array
import os

n = 8
K.<w> = GF(2^n, 'w',repr="log")
KSet = set(K.list())

dual = K.dual_basis()

resultDir = "QAMResults_n=%d"%n
if not os.path.isdir(resultDir):
	os.mkdir(resultDir)

nWorkers = max([4,cpu_count()])

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

print getPolynomFromQAM(H)

def rowSpan(rowVector):
	'''Takes a row vector and outputs the subspace spaned by its elements as a set.'''
	
	N = len(rowVector)
	span = set()
	for ind in combinations[N]:
		newElem = sum([rowVector[i] for i in ind])
		span.add( newElem )
	return span
	
def checkRowRank(row):
	testMatrix = matrix( GF(2), n,0)
	for i in range(row.ncols()):
		testMatrix = testMatrix.augment(matrix(vector(row[0,i])).transpose())
	return testMatrix.rank()
		
		
def checkQAM(mat):
	N = mat.ncols()
	for ind in combinations[N][1:]:		
		testGenerator = matrix(sum(mat[ind,:]))
		if checkRowRank(testGenerator) != N-1:
			print("No QAM because of sum of %s"%str(ind))
			return False
	else:
		print ("QAM")
		return True
		
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
			
####


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
		

def preSearchWorker(solutionList, startMatrix, adaptedDomainFunc, done, quitEv, upperBound):
	print("Process %d starting work."%os.getpid())
	startColumn = matrix(K, startMatrix.nrows(), 1)	
	sol = nextComponent(startColumn, adaptedDomainFunc, startMatrix.nrows()-1)
	
	while not done.is_set():
		solution = sol.next()
		print("Process %d found solution: %s."%(os.getpid(), str(solution)))
		if len(solutionList) < upperBound: solutionList.append(solution)
		if len(solutionList) >= upperBound:
			quitEv.set()
			break
				
def finalSearchWorker(startMatrix, adaptedDomainFunc):
	print("Process %d looks for final QAM."%os.getpid())
	startColumn = matrix(K, n-1, 1)
	nSols = 0	
	N = startMatrix.nrows()
	sol = nextComponent(startColumn, adaptedDomainFunc, N-1)
	while True:
		try:
			solution = sol.next()
			# Save matrix as txt
			B = matrix(K, N + 1, N + 1)
			B[:N,:N] = startMatrix
			B[-1, :N] = solution.transpose()
			B[:N, -1] = solution
			
			with open(resultDir + "/res_%d.txt"%(os.getpid() + nSols), "w") as f:
				f.write(B.str() +"\n" + getPolynomFromQAM(B))
			nSols += 1
		except StopIteration:
			print ("Process %d done, quitting."%os.getpid())
			break
		
### Partition search space
def partitionDomainFunc(originalSetFunc, nProcesses):
	Sf = originalSetFunc
	S0 = list(Sf((0,)))
	S1 = list(Sf((1,)))

	adaptedSetFuncs = []
	if len(S0)/nWorkers >= 1:
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
	
	## Find nWorkers new subQAMs
	A2 = H[:n-2, :n-2]
	S = getOriginalSetFunction(A2)
	adaptedSetFunctions = partitionDomainFunc(S, nWorkers)
	
	m = Manager()
	solutions = m.list()
	workers = []
	quitEv = Event()
	done = Event()
	maxSols = 10*nWorkers
	print("Looking for new subproblems.")
	for partitionFunc in adaptedSetFunctions:
		p = Process(target = preSearchWorker, args = (solutions, A2, partitionFunc, done, quitEv, maxSols))
		p.start()
		workers.append(p)
	quitEv.wait()
	done.set()
	for p in workers:
		p.terminate()
		
	newQAMSubs = []
	A = H[:n-1, :n-1]
	for sol in list(solutions):
		B = copy(A)
		B[-1,:n-2] = sol.transpose()
		B[:n-2,-1] = sol
		print("%s\n"%str(B))
		newQAMSubs.append(B)
		
	#Work on the newly found subproblems 
	workers = []
	nProblems = nWorkers/4		# supposing nWorkers is divisible by 4 such that we can assing 4 cores to each subproblem
	while len(newQAMSubs) > 0:
		for i in range(nProblems):
			sub = choice(newQAMSubs)
			newQAMSubs.remove(sub)
			print("Starting work on \n%s."%str(sub))
			newDomainFuncs = partitionDomainFunc(getOriginalSetFunction(sub), 4)
			for fx in newDomainFuncs:
				p = Process(target = finalSearchWorker, args = (sub, fx) )
				workers.append(p)
				p.start()
			
		for p in workers:
			p.join()
