from constraint import *
from simpleai.search import CspProblem, backtrack
from simpleai.search import MOST_CONSTRAINED_VARIABLE as mcv, HIGHEST_DEGREE_VARIABLE as hdv, LEAST_CONSTRAINING_VALUE as lcv
import numpy as np
import matplotlib.pyplot as plt

def benchmarkRoutine(n):
	K.<w> = GF(2^n, 'w' ,repr="log")
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

	sumIndices = [tuple(ind) for ind in Indices if len(ind) > 1]	
			
	############ Set up CSP using python-constraint ################

	def setUpProblems():
		p = Problem()

		# Add variables and domains
		for i in range(n-1):
			p.addVariable(i, S(i,))

		# Add constraints
		for indices in sumIndices:
			def getConstraintFunc(ind):
				def constraint(*args):
					return sum(args) in S(ind)
				return constraint
			constraintFunc = getConstraintFunc(indices)
			p.addConstraint(constraintFunc, indices)
			
		############ Set up CSP using simpleai ##################
		variables = list(range(n-1))
		domains = [list(S(i,)) for i in range(n-1)]

		constraints = list()
		for indices in sumIndices:
			def getConstraintFunc(ind):
				def constraint(variables, values):
					return sum(values) in S(ind)
				return constraint
			constraints.append((indices,getConstraintFunc(indices)))

		p2 = CspProblem(variables, domains, constraints)

		############ Own recursive iterator #####################

		def getProblemIterator():
			startColumn = matrix(K, n-1, 1)
			
			def getSetFn(oldFn, xi, pos):
				def newSetFn(indices):
					return oldFn(indices).intersection(set([xi + v for v in oldFn(tuple([pos])+indices)]))
				return newSetFn	
					
			def nextComponent(columnTilNow, setFn):
				newCol = copy(columnTilNow)
				pos = n - (newCol.list().count(0) + 1)
				domain = setFn(tuple([pos]))
				if len(domain) > 0:
					for xi in domain:
						newCol[pos] = xi			
						if pos < n-2:
							newSetFn = getSetFn(setFn, xi, pos)
							for col in nextComponent(newCol, newSetFn):
								yield col 
						else:
							yield newCol

			return nextComponent(startColumn, S)
			
		p3it = getProblemIterator()
		
		return p, p2, p3it
		

	####### Benchmarks #########################################

	## Find the first solution
	print("One solution benchmarks.")

	p, p2, p3it = setUpProblems()

	times = dict()
	# python-constraint (pc)
	times["pcBackFC"] = timeit("p.getSolution()", globals={'p':p}, number = 1, repeat = 1, seconds = True)
	p.setSolver(RecursiveBacktrackingSolver())
	times["pcRecBackFC"] =   timeit("p.getSolution()", globals={'p':p}, number = 1 , repeat = 1, seconds = True)

	p.setSolver(BacktrackingSolver(forwardcheck=False))
	times["pcBack"] =    timeit("p.getSolution()", globals={'p':p}, number = 1, repeat = 1, seconds = True)

	p.setSolver(RecursiveBacktrackingSolver(forwardcheck=False))
	times["pcRecBack"] =   timeit("p.getSolution()", globals={'p':p}, number = 1 , repeat = 1, seconds = True)
	
	#p.setSolver(MinConflictsSolver())
	#times["pcMinConf"] = timeit("results['pcMinConf'] = p.getSolution()", number = 1, repeat = 1 , seconds = True)


	# simpleai (sa)
	times["saMCV_Ordered"] =  timeit("backtrack(p2, mcv)",globals={'p2':p2, 'backtrack':backtrack, 'mcv':mcv},  number = 1 , repeat = 1 , seconds = True)
	times["saHDV_Ordered"] =  timeit("backtrack(p2, hdv)", globals={'p2':p2, 'backtrack':backtrack, 'hdv':hdv}, number = 1 , repeat = 1 , seconds = True)
	times["saMCV_LCV"] =  timeit("backtrack(p2, mcv, lcv)", globals={'p2':p2, 'backtrack':backtrack, 'mcv':mcv, 'lcv':lcv}, number = 1, repeat = 1 , seconds = True)
	times["saHDV_LCV"] =  timeit("backtrack(p2, hdv, lcv)", globals={'p2':p2, 'backtrack':backtrack, 'hdv':hdv, 'lcv':lcv}, number = 1, repeat = 1 , seconds = True)
	times["saDefault"] =  timeit("backtrack(p2)",globals={'p2':p2, 'backtrack':backtrack}, number = 1 , repeat = 1, seconds = True)

	# my implementation
	times["own"]=   timeit("p3it.next()", globals={'p3it':p3it},number = 1, repeat = 1, seconds = True)

	## Find 5 Solutions
	print("Multiple solution benchmarks.")
	timesMult = dict()
	p, p2, p3it = setUpProblems() 	# reset search trees
	pit = p.getSolutionIter()

	timesMult["pcFC"] =   timeit("[pit.next() for i in [1,2,3,4,5]]", globals = {'Integer':Integer,'pit':pit}, number = 1, repeat = 1, seconds=True)
	p.setSolver(BacktrackingSolver(forwardcheck = False))
	pit = p.getSolutionIter()
	timesMult["pc"] =   timeit("[pit.next() for i in [1,2,3,4,5]]", globals = {'Integer':Integer,'pit':pit}, number = 1, repeat = 1, seconds = True)
	timesMult["own"] =   timeit("[p3it.next() for i in [1,2,3,4,5]]", globals = {'Integer':Integer,'p3it':p3it}, number = 1 ,repeat = 1, seconds = True)

	#print(rankingMult)
	
	return times, timesMult
	
# Do benchmarks for different n and plot everything

for i in range(7,8):
	print ("Benchmarking dimension %d."%i)
	dataSingle, dataMulti = benchmarkRoutine(i)
	
	###### Single solution
	vals = dataSingle.values()
	# scaling 
	fastest = min(vals)
	factor = 1
	power = 0
	while fastest < 1:
		factor = factor * 10
		fastest = fastest * 10
		power += 1
		
	vals = [factor*v for v in vals]
	ylabel = "Time in sec." if factor == 1 else "Time in sec.*E-%d"%power
	tickLabels = dataSingle.keys()
	ind = np.arange(len(tickLabels))
	
	fig, ax = plt.subplots()
	bars = ax.bar(ind, vals, width = 0.2)
	ax.set_ylabel(ylabel)
	ax.set_xlabel("Solver")
	ax.set_title("'Last Column Search', 1 Sol. in %d dimensions."%i)
	ax.set_ylim(top=1.15 * max(vals))
	plt.xticks(ind, tickLabels, rotation="vertical")
	plt.subplots_adjust(bottom=0.3)
	
	for rect in bars:
		height = rect.get_height()
		ax.text(rect.get_x() + rect.get_width(), 1.05*height,'%.2f' % height, ha='center', va='bottom')
	plt.tight_layout()
	plt.savefig("barChartSingleSol%d.png"%i)
	plt.close()
	
	######## 5 Solutions
	vals = dataMulti.values()
	# scaling 
	fastest = min(vals)
	factor = 1
	power = 0
	while fastest < 1:
		factor = factor * 10
		fastest = fastest * 10
		power += 1
		
	vals = [factor*v for v in vals]
	ylabel = "Time in sec." if factor == 1 else "Time in sec.*E-%d"%power
	tickLabels = dataMulti.keys()
	ind = np.arange(len(tickLabels))
	
	fig, ax = plt.subplots()
	bars = ax.bar(ind, vals, width = 0.4)
	ax.set_ylabel(ylabel)
	ax.set_xlabel("Solver")
	ax.set_title("'Last Column Search', 5 Sol. in %d dimensions."%i)
	ax.set_ylim(top=1.15 * max(vals))
	plt.xticks(ind, tickLabels, rotation="vertical")
	plt.subplots_adjust(bottom=0.3)
	
	for rect in bars:
		height = rect.get_height()
		ax.text(rect.get_x() + rect.get_width(), 1.05*height,'%.2f' % height, ha='center', va='bottom')
	plt.tight_layout()
	plt.savefig("barChartMultiSol%d.png"%i)
	plt.close()
