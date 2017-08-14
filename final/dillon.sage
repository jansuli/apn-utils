n = 3

# Set up fields and vector spaces
F2 = GF(2)
V = VectorSpace(F2, n)
VBasis = V.basis()					# column options for L in row echelon
K.<w> = GF(2^n, 'w', repr="log") 				# w as primitive Element
KList = K.list()					# for comfortable indexing
CodeSpace = VectorSpace(F2, 2^n-1)
R.<z> = PolynomialRing(K, 'z')		# needed for Lagrange-Interpolation, thus for actually printing polynomials (in z)

def lagrangeInterpolation(func):
    poly = R(0)
    for i in range(2^n):
        ai = KList[i]
        bi = func(ai)
        prod = R(1)
        for k in range(2^n):
            if k != i:
                ak = KList[k]
                prod = prod* ((z-ak)/(ai-ak))
        poly = poly + bi*prod
    return poly

# Test kim function 
def f(x):
	return x^3 + x^10 + w* x^24
	
# Build generator matrix H of dual Code
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

Hext = matrix(F2, [0 for i in range(2*n)]).transpose().augment( H ) 	# H with added zero column
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
			zeta[trailingZeros].append(c)#matrix(c).transpose())		# store as transposed matrix already to speed up search later on
		
## Find all simplex sub codes 	
solutions = list()
print("Looking for simplex-sub-codes.")
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
		testVectors += [c[:j] for c in zeta[i]]
	
	for option in options:
		testMatrix = mat.augment(option)
		for c in testVectors:
			if testMatrix * c == 0:
				break
		else:
			if j <= 2*n - 1:
				nextColumn(testMatrix)
			else:
				print("Solution:\n%s\n"%testMatrix.str())
				solutions.append(testMatrix)
			
startMatrix = matrix(F2, n, 0)
nextColumn(startMatrix)

## Search disjunkt spaces
pairs = list()
print("Found %d simplex-sub-codes. Now checking %d combinations for disjoint pairs."%(len(solutions), binomial(len(solutions),2))) 
for pair in Combinations(solutions, 2):
	L1 = pair[0]
	L2 = pair[1]
	C1 = CodeSpace.subspace( (L1*H).rows() )
	C2 = CodeSpace.subspace( (L2*H).rows() )
	
	if C1.intersection(C2).dimension() == 0:
		pairs.append(pair)
		
## If we found at least one disjoint pair, print the first apn-permutation
if len(pairs) > 0:
	print("We have %d disjoint pairs. The first one leads to the following APN-permutation:"%(len(pairs)))
	L1, L2 = pairs[0]
	H1 = L1 * Hext
	H2 = L2 * Hext
	
	def pi2(x):
		pos = KList.index(x)
		resultVector = H2.columns()[pos]
		return K(resultVector)
		
	def pi1Inv(x):
		vec = vector(x)
		pos = H1.columns().index(vec)
		return KList[pos]
		
	def g(x):
		return pi2( pi1Inv(x) )
		
	gPoly = lagrangeInterpolation(g)
	print("g(x) = "+str(gPoly))
	
