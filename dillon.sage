try:
	np = numpy
	pk = pickle
	tqdm = tqdm
	sleep = slepp
except NameError:
	import numpy as np
	import pickle as pk
	import sys
	from tqdm import tqdm
	from time import sleep

#from __future__ import print_function

m=6
K.<a> = GF(2^m, 'a')
print("The minimal polynomial of 'a' is:")
print(a.minimal_polynomial())
def f(x):
  return x^3
  
def g(x):
	return x^3 + x^10 + a * x^24

def generatorGF(K, f):
  '''Returns the 2x2^m-1 generator matrix H_f of Code C_f as defined by Dillon.
  Takes a finite field K and an APN-function f.'''

  H = []
  a = K.primitive_element()
  m = K.degree()
  
  for i in range(0, 2^m-1):
    H.append([a^i, f(a^i)])
  return matrix(K, H).transpose()

def GFtoBinMatrix(M, dim):
  # iterate over the columns
  A = []
  N = M.ncols()
  for i in range(0, N):
    col = M[:, i]
    tmpCol = []
    for elem in col.list():
      coeff = vector(elem).list()
      tmpCol = tmpCol + coeff
    A.append(tmpCol)
  return matrix(GF(2),A).transpose()

def evalsZero(mat, listOfVectors):
	for vec in listOfVectors:
		#print("Teste den Vector {}".format(vec))
		if mat * vec == 0:
			if mat[:,0] == 0: tqdm.write(mat.str() + "*"+ str(vec.list()) +"evals to zero")
			MAT = str(mat*vec)
			if mat[:,0] == 0: tqdm.write(MAT +"\n\n") 
			return True

def altSimplexMatrices(ZDict, K):
	m = K.degree()
	m2 = m*2
	V = K.vector_space()
	v = V.list()
	v_sorted,vClasses = sortedSumSet(matrix(v).transpose(),verbosity=True, sumUp = False)
	start = vClasses[m] + vClasses[m-1] # due to reduced row echelons
	baseMatrixV = matrix.identity(GF(2),m)
	
	testVectorLookUp = dict()
	
	for j in range(1,m2):
		testVecs = list()
		for i in range(m2 - j-1, m2+1):
			print("Including class %d into testVecs."%i)
			testVecs = testVecs + ZDict[i]
		testVecs = [vec[:j+1] for vec in testVecs]
		testVectorLookUp[j] = testVecs
		print "\n\n\n"
		
	
	def nextCol(matrixTilNow):
		nCols = matrixTilNow.ncols()
		rank = matrixTilNow.rank()
		lastCol = matrixTilNow[:,-1].list()
		lastClass=lastCol[::-1].index(1) if 1 in lastCol else m
		requiredCols = m2-nCols
		lookUp = list() 
				
		# Vectors in v, that allow row reduced echelon form:
		if lastClass != 0 and (m-rank) < requiredCols:
			for i in range(lastClass, m+1):
				lookUp = lookUp + vClasses[i]
			lookUp.append(vector(baseMatrixV[-lastClass]))
		elif lastClass != 0:
			lookUp = [  vector(baseMatrixV[-lastClass])]
		else:
			lookUp = v
		#tqdm.write(matrixTilNow.str() + "\n")
		
		if(len(lookUp) > 0):
			testVecs = testVectorLookUp[nCols]
			
			for candidate in lookUp:#, postfix={"sols: ":len(sols)}):
				candMatrix = matrixTilNow.augment(candidate)
				print ("While havin %d solutions investigating \n%s ...\n\n"%(len(sols),candMatrix.str()))
				kernel = candMatrix.right_kernel()
				for testVec in testVecs:
					if testVec in kernel:
						verworfen.append(candMatrix)
						break
				else:
					if requiredCols != 1:
						nextCol(candMatrix)
					else:
						sols.append(candMatrix)
						print(len(sols))
	sols=list()
	verworfen = list()
	start.reverse()
	for vec in start:
		mat = matrix([vec]).transpose()
		nextCol(mat)
	
	return sols 
	
def altSimplexMatricesD(ZDict, K):
	m = K.degree()
	m2 = m*2
	V = K.vector_space()
	v = V.list()
	v_sorted,vClasses = sortedSumSet(matrix(v).transpose(),verbosity=True, sumUp = False)
	start = vClasses[m] + vClasses[m-1] # due to reduced row echelons
	baseMatrixV = matrix.identity(GF(2),m)
	
	testVectorLookUp = dict()
	
	for j in range(1,m2):
		testVecs = list()
		for i in range(m2 - j-1, m2+1):
			print("Including class %d into testVecs."%i)
			testVecs = testVecs + ZDict[i]
		testVecs = [vec[:j+1] for vec in testVecs]
		testVectorLookUp[j] = testVecs
		print "\n\n\n"
		
	
	def nextCol(matrixTilNow):
		nCols = matrixTilNow.ncols()
		rank = matrixTilNow.rank()
		lastCol = matrixTilNow[:,-1].list()
		requiredCols = m2-nCols
		lookUp = list() 
				
		# Vectors in v, that allow row reduced echelon form:
		if rank != m and (m-rank) < requiredCols:
			for i in range(m-rank, m+1):
				lookUp = lookUp + vClasses[i]
			lookUp.append(vector(baseMatrixV[rank-m]))
		elif rank != m :
			lookUp = [  vector(baseMatrixV[rank-m])]
		else:
			lookUp = v
		lookUp.reverse()
		#tqdm.write(matrixTilNow.str() + "\n")
		
		if(len(lookUp) > 0):
			testVecs = testVectorLookUp[nCols]
			for candidate in lookUp:#, postfix={"sols: ":len(sols)}):
				candMatrix = matrixTilNow.augment(candidate)
				print ("While havin %d solutions investigating \n%s...\n\n"%(len(sols),candMatrix.str()))
				kernel = candMatrix.right_kernel()
				for testVec in testVecs:
					if testVec in kernel:
						break
				else:
					if requiredCols != 1:
						nextCol(candMatrix)
					else:
						sols.append(candMatrix)
						print(len(sols))
	sols=list()
	verworfen = list()
	start.reverse()
	for vec in start:
		mat = matrix([vec]).transpose()
		nextCol(mat)
	
	return sols,verworfen 
	
def findSimplexMatrices(ZDict, K, verbosity=True, save=True, reduced_echelon = True):
	m = K.degree()
	m2 = m*2
	V = K.vector_space()
	v = [e.list() for e in V]
	sols=list()
	v_sorted,vClasses = sortedSumSet(matrix(v).transpose(),verbosity=True, sumUp = False)
	for i in range(0, m+1):
		vClasses[i] = [vec.list() for vec in vClasses[i]]
	
	progressBarLabels = dict()
	for i in range(0, m2):
		progressBarLabels[i] = 0
	
	def findNextCol(cols):
		if verbosity == 2: tqdm.write("Start vector was:")
		if verbosity == 2: tqdm.write(vector(cols[0]).str())
		nCols = len(cols)
		classN = cols[-1][::-1].index(1) if 1 in cols[-1] else m2		# class of last constructed column
		reqN = m2 - nCols												# columns left to build. at the moment constructing col no. len(cols) + 1
		
		def generateLookUp(vecs,colsTilNow):
			rank = matrix(colsTilNow).rank()
			lastCol = colsTilNow[-1]
			lastClass = lastCol[::-1].index(1) if 1 in lastCol else m
			if verbosity == 2: tqdm.write("last class: %d"%lastClass)
			### pre-filtering: A we need m steps only vectors of classes lastClass to m and the e of lastClass-1are allowed:
			lookUp=list()
			if lastClass != 0 and (m - rank) < reqN:
				for i in range(lastClass if reduced_echelon else lastClass-1, m+1):
					lookUp = lookUp + vClasses[i]
				
				if reduced_echelon:
					e = [0 for i in range(0,m)]		## strenge zeilenstufen
					e[m-lastClass] = 1
					lookUp.append(e)
			elif lastClass != 0:
				e = [0 for i in range(0,m)]		## strenge zeilenstufen
				e[m-lastClass] = 1
				lookUp =[e]
			else:
				lookUp = vecs
			
			return lookUp
				
		# check for row echelon and rank
		lookUp = generateLookUp(v, cols)
		progressBarLabels[nCols] += 1
	
		if verbosity == 2: tqdm.write("We have %d possible vectors to choose from and test."%(len(lookUp)))
		
		
		if(len(lookUp)>0):
			# Testing against sorted sum set given by ZDict
			testVecs = ZDict[m2-nCols-1]
			#for i in range(m2,m2-nCols-2, -1):
			#	testVecs = testVecs + ZDict[i]
			testVecs = [vec[:nCols+1] for vec in testVecs]
			for vec in tqdm(lookUp, postfix = {'#sols': len(sols)}, desc="level %d.%d"%(nCols, progressBarLabels[nCols])):
				newCols = list(cols)
				newCols.append(vec)
				mat = matrix(newCols).transpose()
				
				if not evalsZero(mat, testVecs):
					#print("possible new col :D Until now we found %d complete solutions"%len(sols))
					if(reqN > 1 ):
						findNextCol(newCols)
					elif(reqN == 1):
						if verbosity == 2: tqdm.write("SUCCESS, now we got #sols: %d"%(len(sols)))
						if verbosity == 1: tqdm.write(mat.str() + "\n \n")
						sols.append(mat)
						if save :
							with open("sols", "w+") as f:
								pk.dump(sols, f)
						
						#return True # remove later on
				elif(verbosity == 2):
					tqdm.write("The following matrix evaluates to zero:")
					tqdm.write(mat.str())
					tqdm.write("------------------------------------------------\n")
		else:
			tqdm.write("Not enough choices left.")
		
	#test run with one start vector choosen at random
	#vecs.remove(vecs[0])
	start = [[0 for i in range(0,m)]]
	e = [0 for i in range(0,m)]
	e[0] = 1
	start.append(e)
	for vec in start:
		mat = [vec]
		findNextCol(mat)
	#mat = [vecs[1]]
	#findNextCol(mat) 
	return sols
	
def checkDisjoint(matrixList, H):
	V = VectorSpace(GF(2), H.ncols())
	disjointPaires = []
	for mat1 in tqdm(matrixList):
		for mat2 in matrixList:
			if mat1 != mat2:
				Mat1 = mat1*H
				Mat2 = mat2*H
				
				V1 = V.subspace(Mat1.rows())
				V2 = V.subspace(Mat2.rows())
				
				if V1.intersection(V2).dimension() == 0:
					tqdm.write("Found disjoint subspaces!!! \n")
					disjointPaires.append((matrixList.index(mat1), matrixList.index(mat2),[Mat1, Mat2]))
				else:
					tqdm.write("Subspaces not disjoint")
	return disjointPaires
	
def sortedSumSet(H, verbosity = False, sumUp = True):
	M = H.ncols()
	if sumUp:
		vecSums = list ()		# save new vectors
		for i in range(0,M):
			x = H[:, i]			# first column vector
			for j in range(0,M):
				y = H[:, j]		# second cocumn vector
				if (x != y):
					z = (x+y).list()
					if z not in vecSums:		# if results is not yet in list, save it there
						vecSums.append( z )
	else:
		vecSums =H.columns()
	if verbosity: print("Our set Σ has cardinality %d."%len(vecSums)) 
		
	def sortByRow(B,rowIndex):
		''' Takes numpy array. '''
		currentClass = nRows - rowIndex - 1
		B = np.array(B) 	# copying
		[l,k] = B.shape
		if verbosity: print("Entered function to look for vectors of class %d in a %d×%d Matrix."%(currentClass, l,k))
		if k != 0:			# do we have colums left to sort?
			rowB = B[rowIndex, :].tolist()
			sorting = np.argsort(rowB)
			B[:,:] = B[:, sorting]
			rowB = B[rowIndex, :].tolist()
			#print("Sorted row for rowIndex %d is: "%(rowIndex), rowB)
			if 1 in rowB:
				onePos = rowB.index(1)
				if verbosity: print("There are vectors of class %d in Σ."%(currentClass))
				possibleClasses.remove(currentClass)
				subR = B[:, onePos:]
				subL = B[:, :onePos]
				vecSumClassDict[currentClass] = matrix(GF(2),subR).columns()
				if verbosity: print("These are:")
				
				if verbosity: print(subR)
				
				if rowIndex > 0:
					subLS = sortByRow(subL, rowIndex - 1)
					return np.hstack((subLS, subR))
				else:
					
					if subL.shape[1] != 0:
						print("Class m included")
						possibleClasses.remove(nRows)
						vecSumClassDict[nRows] = matrix(GF(2),subL).columns()
					return np.hstack((subL, subR))
			else:
				if verbosity: print("There are no vectors of class %d in Σ. Checking for next..."%(currentClass))
				if rowIndex > 0:
					return sortByRow(B, rowIndex-1)
				else:
					return B
		else:
			if verbosity: print("There is nothing left to sort. Returning.")
			return B
	# Start recursion by sorting A by last Row:
	vecSumClassDict = dict()		# will be filled partially by sortByRow
	A = np.array(vecSums).transpose()
	nRows = A.shape[0]
	possibleClasses = [i for i in range(0, nRows + 1)]	# will be modified by sortByRow
	A_sorted = sortByRow(A, nRows-1)
	
	# Assign empty list to empty classes:
	for classN in possibleClasses:
		vecSumClassDict[classN] = []
	return A_sorted, vecSumClassDict

print("Generator Matrix H over finite Field GF(2^m)")
H = generatorGF(K,f)
print("Binary Generator h.")
h = GFtoBinMatrix(H, m)


print("Creating sorted sum set...")
#z, zeta = sortedSumSet(h.augment(vector([0 for i in range(0,2*m)])), verbosity =True)

#sols = findSimplexMatrices(zeta, K, verbosity=1) 
#sols = altSimplexMatrices(zeta, K)
#dis = checkDisjoint(sols, h)

def f1(x):
	b = a^(-2)
	return b^38*x^48+b^33*x^40 + b^28*x^34 + b^25*x^33 + b^43 * x ^ 32 + b^5* x^24 + b^42 *x^20 + x^17 + b^2*x^16 + b^4* x^12 + b^7 *x^10 + b^58*x^8 + b^59 *x^6 + b^5 *x^5 + b^36 * x^4 + b^47 *x^3 + b^30 *x^2 +b^9 *x
	
#F1 = generatorGF(K, f1)[1,:]
#F1b = GFtoBinMatrix(F1, m)

#F1L = h.solve_left(F1b).echelon_form()

	
#debugSimplex(zeta,K,F1L)

#sols, verw = altSimplexMatricesD(zeta,K)
#checkDisjoint(sols)
