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
		
		
