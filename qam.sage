from constraint import *

n = 7
K.<w> = GF(2^n, 'w' ,repr="log")

KBasis = [w^i for i in range(n)]
DualBasis = K.dual_basis()
CanBasis = VectorSpace(GF(2), n-1).basis()

comb = Combinations(n-1).list()
combN = Combinations(n).list()[1:]	# nontrivial row combinations

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
	print testMatrix.str()
	return testMatrix.rank()

def checkQAM(mat):
	for ind in combN:		
		testGenerator = matrix(sum(mat[ind,:]))
		print testGenerator
		if checkRowRank(testGenerator) != n-1:
			print("no qam")
			return False
	else:
		print ("QAM")
		return True

def rowSpan(row, matrix = False):
	ret = set()
	if matrix:
		print "Matrix"
		for ind in comb:
			ret.add( sum( row[0, ind].list() ))
		return ret
	for ind in comb:
		ret.add( sum([row[i] for i in ind]) )
	return ret	
	
def f(x):
	return x^3
	
Cf = matrix(K,n,n)
Cf [0,1] = Cf[1,0] = 1

H = M.transpose()*Cf*M
print H.str() + "\n"

A = H[:n-1, :n-1]

KSet = set(K.list())

setFunctions = list()
def Sf(indexVec):
	#print("First set function")
	return KSet.difference( rowSpan(indexVec*A))
setFunctions.append(Sf)


## Demonstrate how to work with set functions:
#i = 0
#S1 = setFunctions[i](CanBasis[i])

##checkQAM(H)

#options = dict()
#options[i] = list(S1)
#options[n-2] = []
#i = 0
#while options[n-2] == []:
	#for xi in options[i]:
		#i = 0
		#used = [K(0)]
		#H = copy(Hor)
		#while i < n-2:
			##print("Dealing with %s"%str(xi))
			#used.append(xi)
			#H[i, n-1] = H[n-1, i] = xi
			##print("%s\n"%H.str())

			#def defineSetFunc(i):
				#def Sf(indexVec):
					##print("Generating set for i = %d"%(i+1))
					#return setFunctions[i](indexVec).intersection(set([s-xi for s in setFunctions[i](CanBasis[i] + indexVec)]))
				#return Sf
			#sf = defineSetFunc(i)
			#setFunctions.append(sf)
			#S2 = sf(CanBasis[i])
			#lS2 = list(S2)
			#options[i + 1] = lS2
			##print(len(S2),i)
			#newOptions = list(S2.difference(set(used)))
			#if len( newOptions ) > 0:
				#xi = choice(newOptions)
			#else:
				#print("no choices at %d"%i)
				#break
			#i += 1
		##print "done"
		#if i == n-2:
			#print "SUCCESS"
			#H[n-1, n-2] = H[n-2, n-1] = xi
			#break

# Treat as CSP

domains = dict()
for i in range(A.nrows()):
	domains[i] = list(Sf(CanBasis[i]))
	
p = Problem()
for i in range(A.nrows()):
	p.addVariable(i, domains[i])
	
sumIndices = [ind for ind in comb if len(ind) > 1]

def getConstraint(avoid):
	def SumNotInSet(*args):
		return sum(args) not in avoid
	return SumNotInSet 
	
for ind in sumIndices:
	indexVec = matrix(GF(2), 1, n-1)
	indexVec[0, ind] = 1
	#print indexVec
	avoidSet = list(rowSpan((indexVec*A).list()))
	constraintFunction = getConstraint(avoidSet)
	p.addConstraint(constraintFunction,ind)
	print ("Sum of rows %s constraint."%str(ind))

def applySolution(mat, solDict):
	for i in solDict:
		mat[n-1, i] = mat[i, n-1] = solDict[i]
		
it = p.getSolutionIter()
while True:
	try:
		sol = it.next()
		Hmod = copy(H)
		applySolution(Hmod,sol)
		print ("%s\nis a solution."%Hmod.str())
	except StopIteration:
		print ("No more solutions.")
	

	
