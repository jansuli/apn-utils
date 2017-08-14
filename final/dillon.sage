n = 6

# Set up fields and vector spaces
F2 = GF(2)
V = VectorSpace(F2, n)
VBasis = V.basis()
K.<w> = GF(2^n, 'w') 	# w as primitive Element
CodeSpace = VectorSpace(F2, 2^n-1)

# Test kim function 
def f(x):
	return x^3 + x^10 + w* x^24
	
# Build generator matrix of dual Code
H = matrix(F2, 2*n, 0)
for x in K.list()[1:]:
	top = vector(x).list()
	bottom = vector(f(x)).list()
	column = matrix(F2, top + bottom).transpose()
	H = H.augment(column)
	
# Possible columns of L, sorted by number of trailing zeros
# (needed to build L in row reduced echelon form)
columnOptions = dict()
for i in range(n + 1):
	columnOptions[i] = list()
	
for vec in V:
	vecList = vec.list()
	vecList.reverse()
	if 1 in vec:
		trailingZeros = vecList.index(1)	# index of first 1 in reversed vector corresponds to number of trailing zeros, index=0 => no trailing zeros etc.
	else:
		trailingZeros = n
	columnOptions[trailingZeros].append(vec)
	
# Build set of column sums Σ
# again sorted by number of trailing zeros to get the test equations
zeta = dict()
for i in range(2*n + 1):
	zeta[i] = list()

Hext = H.augment( matrix(F2, [0 for i in range(2*n)]).transpose() )
for x in Hext.columns():
	for y in Hext.columns():
		if x[:n] != y[:n]:
			c = x + y
			cList = c.list()
			cList.reverse()
			
			if 1 in cList:
				trailingZeros = cList.index(1)
			else:
				trailingZeros = 2*n
			zeta[trailingZeros].append(matrix(c).transpose())
		
## Find all simplex sub codes 	
solutions = list()

def nextColumn(mat):
	j = mat.ncols()	+ 1		# column to construct, j as explained in thesis
	r = mat.rank()
	# to get row reduced echelon form we may take any vector with number of trailing zeros greater than n - r and the (r+1)st basis vector
	options = [] if r == n else [ VBasis[r] ]
	for i in range(n-r, n+1):
		options += columnOptions[i]
	
	# test against all vectors in Σ with number of trailing zeros at least 2*n - j
	testVectors = []
	for i in range(2*n-j, 2*n+1):
		testVectors += [c[:j, 0] for c in zeta[i]]
	
	for option in options:
		testMatrix = mat.augment(option)
		print("Testing \n%s"%(testMatrix.str()))
		for c in testVectors:
			if testMatrix * c == 0:
				break
		else:
			if j <= 2*n - 1:
				nextColumn(testMatrix)
			else:
				print("Solution:\n%s"%testMatrix.str())
				solutions.append(testMatrix)
			
startMatrix = matrix(F2, n, 0)
nextColumn(startMatrix)

## Search disjunkt spaces
pairs = list()
for pair in Combinations(solutions, 2):
	L1 = pair[0]
	L2 = pair[1]
	C1 = CodeSpace.subspace( (L1*H).rows() )
	C2 = CodeSpace.subspace( (L2*H).rows() )
	
	if C1.intersection(C2).dimension() == 0:
		pairs.append(pair)
