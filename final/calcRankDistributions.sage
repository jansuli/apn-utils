from multiprocessing import Process, Manager, cpu_count
from Queue import Empty
import pickle
from os import getpid

from funcs import functions 	# just a big function string

class APN():
	'''Class representing quadratic APN functions F^n -> F^n.'''
	
	def __init__(self, field, expression, primName = 'a'):
		self.field = field			
		self.vecSpace = field.vector_space()
		self.expression = expression
		
		self.primName = primName	# name of primitive element used in expression string
		self.primElem = field.primitive_element()
	
	def func(self,x):	
		'''Evaluate function. Takes field element as input and returns field element.'''
		
		return sage_eval(self.expression,locals = {'x':x, self.primName: self.primElem} )
		
	def vecFunc(self, vec):
		'''Evaluates function taking a vector as input and outputs vector, too.'''
		
		fieldElem = self.field(vec)
		return vector(self.func(fieldElem))
		
	def componentFunctions(self):
		if not hasattr(self, "componentFuncs"):
			componentFuncs = dict()
			def getComponentFunction(b):
				def Fb(v):
					res = b.dot_product(self.vecFunc(v))
					return res
				return Fb
				
			# Iterate over all codomain vectors to generate compontent functions.
			for b in self.vecSpace:
				b.set_immutable()	# set immutable so it can be used as a key
				componentFuncs[b] = getComponentFunction(b)
			self.componentFuncs = componentFuncs
		
		return self.componentFuncs
		
def alternatingBilinearForm(func, vecSpace):
	'''Maps boolean function 'func':F^n -> F (we will later use the component functions) to alternating bilinear form.'''
	
	def Bf(x,y):
		return func(x+y) + func(x) + func(y) + func(vecSpace.zero())
	
	return Bf
	
def alternatingMatrix(form, basis):
	'''Takes an alternating bilinearform F^n × F^n -> F and returns its binary n×n matrix.'''
	
	mat = matrix(GF(2), dim)
	for i in range(1, dim):
		for j in range(0,i):
			mat[i,j] = form(basis[i], basis[j])
	mat = mat + mat.transpose()
	return mat

nWorkers = cpu_count()
	
dim = 8
K.<w> = GF(2^dim, 'w')
V = K.vector_space()
VBasis = V.basis()

def calcRankDist(jobs, resultQueue):
	print("Process %d started work."%getpid())
	while True:
		try:
			apnString = jobs.get_nowait()
			f = APN(K, apnString, 'g')
			
			componentFunctions = f.componentFunctions()
			
			print("Iterating over component functions...")
			rankCount = dict()
			for i in range(5):
				rankCount[i] = 0	# possible ranks: 0,2,4,6,8 -> indexing 0,1,2,3,4
				
			for b,Fb in componentFunctions.iteritems():
				Bf = alternatingBilinearForm(Fb, V)
				B = alternatingMatrix(Bf, VBasis)
				r = B.rank()
				rankCount[r/2] +=1
				
			rankDistribution = tuple([rankCount[r] for r in range(5)])
			print("For %s… the rank distribution is %s."%(apnString[:15], str(rankDistribution)))
			
			resultQueue.put( (apnString, rankDistribution) )
		except SyntaxError:
			rankQueue.put( (apnString, "Threw an error...") )
		except Empty:
			print("No more functions to explore. Worker %d quitting."%getpid())
			break

if __name__ == "__main__":
	manager = Manager()
	rankDists = manager.Queue()
	jobs = manager.Queue()
	
	functionStrings = functions.replace("\n","").split(',')[:2]
	for funcString in functionStrings:
		jobs.put(funcString)
	
	workers = []
	print("Investigating %d function strings."%jobs.qsize())
	
	for i in range(nWorkers):
		p = Process(target=calcRankDist, args=(jobs, rankDists))
		p.start()
		workers.append(p)
	for p in workers:
		p.join() 	# wait for each process to finish
		
	print("\nDone with work. Now saving the results.")
	print("This is done via pickling a dict with distribution-tuples as keys and lists of corresponding functions as values.")
	
	results = dict()
	while not rankDists.empty():
		apnStr, rankDistribution = rankDists.get()
		if not rankDistribution in results:
			results[rankDistribution] = [apnStr]
		else:
			results[rankDistribution].append(apnStr)
			
	with open("manyRankDists.data", "w") as f:
		pickle.dump(results, f)
		
