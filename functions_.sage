def GFtoBinMatrix(M):
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
  
def checkApnMatrix(mat):
	''' Taking a the parity check matrix 'mat' of a Code this function returns True if the Code has minimum distance d, else False.
	
	It does so by checking each tuple of 4 different columns for linear independce (via matrix rank) and looking for a 5-tuple of linear dependent columns.
	This process can get quite tedious as for the first step alone we have to check binomial{mat.ncols()}{4} combinations.'''
	d = 5
	n = mat.ncols()
	indices = list()
	for a in tqdm(range(n)):
		for b in range(n):
			if a!= b:
				for c in range(n):
					if a!=b!=c:
						for e in range(n):
							if (a != b != c != e): indices.append( [a,b,c,e] )
	print len(indices)
	D = d-1
	print("Looking for indipendent %d columns."%(D))
	for ind in tqdm(indices):
		r = mat[:, ind].rank()
		
		if r != D:
			print("failed for indices %s."%str(ind))
			return False
	else:
		print("Every possible combination of 4 different columns is lineary independent.")
		
	indices = Combinations(range(n), d).list()
	for ind in tqdm(indices):
		M = mat[:, ind]
		
		if M.rank() != 5:
			print("Found a combination of d lin. dependent columns:\n" + M.str())
			return True
	else:
		return False
  
def alternatingBilinearForm(func, vecSpace):
	'''Maps boolean function 'func':F^m -> F (we will later use the component functions) to alternating bilinear form.'''
	
	def Bf(x,y):
		return func(x+y) + func(x) + func(y) + func(vecSpace.zero())
	return Bf
	
def alternatingMatrix(form, basis):
	size = basis[0].parent().dimension()
	mat = matrix(GF(2), size)
	for i in range(1,size):
		for j in range(0,i):
			mat[i,j] = form(basis[i], basis[j])
	mat = mat + mat.transpose()
	return mat
	
def rankDistCleanUp(rankDist):
	m = len(rankDist) - 1
	k = floor(m/2)
	newRankDist = dict()
	for i in range(0, k + 1):
		newRankDist[i] = rankDist[2*i]
	return newRankDist
			

def codeDistanceFromRanks(rankDist,m):
	'''Takes a CLEANED UP rank distribution dict (and dimension m) and returns weight distribution of corresponding RM code CF.'''
	weightDist = dict()
	k = len(rankDist) - 1
	for i in range(0, 2^m + 1):
		if i%2 == 0:
			if i == 0 or i == 2^m:
				weightDist[i] = 1
			elif i == 2^(m-1):
				weight = 2^(2*m+1)
				for h in range(0, k+1):
					weight -= rankDist[h] * 2^(2*h)
				weightDist[i] = weight
			else:
				for h in range(1,k+1):
					if i == 2^(m-1) + 2^(m-h-k) or i == 2^(m-1) - 2^(m-h-k):
						weightDist[i] = rankDist[h] * 2^(2*h)
		else:
			pass#weightDist[i] = 0
	return weightDist
	
def WalshFromRanks(rankDist,m):
	walsh = list()
	k = len(rankDist) - 1
	for i in range(0,2^m +1):
		if i%2==0:
			if i == 2^m:
				walsh.append(i)
			elif i == 0:
				exp = 2^(2*m)
				for h in range(0,k+1):
					exp -= rankDist[h]*2^(2*h)
				walsh.append(i^exp)
			else:
				for h in range(1,k+1):
					if i == 2^(m-h):
						walsh.append(i^(rankDist[h]*2^(2*h)))
		else:
			walsh.append(0)
	return walsh
		
def checkApnByDefinitionDict(funcTable, K):
	k = list(K.list())
	#k.remove(K(0))
	
	def newFunc(a):
		def f(x):
			return funcTable[x+a] - funcTable[x]
		return f

	for a in tqdm(k):
		f = newFunc(a)
		for b in K:
			sols = 0
			for x in K:
				if f(x) == b:
					sols += 1
			if sols > 2:
				#tqdm.write( a,b, sols)
				return False
	return True
			
def checkApnByDefinition(func, K):
	
	def returnDerivative(a):
		def Df(x):
			return func(x+a) - func(x)
		return Df

	for a in tqdm(K, desc="a-iterator"):
		if a != 0:
			f = returnDerivative(a)
			for b in K:
				sols = 0
				for x in K:
					if f(x) == b:
						sols += 1
				if sols > 2:
					print( a,b, sols)
					return False
	return True
	
def uniKim( (a,b,c) ):
	def kim(x):
		return sage_eval("x^3 + a*x^32 + b*x^65 + c*x^96", locals={'a':a,'b':b,'c':c, 'x':x})
	return kim
	
M = matrix(K, m)
for i in range(10):
	for j in range(10):
		M[i,j] = (w^j)^(2^i)
	
def thouroghCheckQAM( (a,b,c) ):
	
	C = matrix(K, m)

	C[5,1] = C[1,5] = a
	C[6,0] = C[0,6] = b
	C[6,5] = C[5,6] = c

	H = M.transpose() * C * M
	
	cols = H.columns()
	B = V.subspace(cols)
	for vec in B:
		if vec != 0:
			rCheck = matrix(vector(vec[0]))
			for elem in vec[1:]:
				rCheck = rCheck.stack(vector(elem))
			r = rCheck.rank()
			if r != m-1:
				tqdm.write(str(r))
				return False
	else:
		return True
		
