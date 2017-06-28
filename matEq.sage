from tqdm import tqdm
k = 3
n = 2^k-1
m = binomial(n,4)*15
G.<y> = GF(2^(2*k), 'y')
K.<w> = GF(2^k, 'w')
V = VectorSpace(G, m)
#K = G.subfields(m, 'w')[0][0]

def f(x):
	return x^3
	
basis = [y^i for i in range(2*k)]
def kVecTok2field( vec, offset=0):
	elem = G(0)
	for i in range(len(vec)):
		if vec[i] != 0:
			elem += basis[i+offset]
	print("Transformed vector %s into field element %s (with offset %d)."%(str(vec), str(elem), offset))		
	return elem

print("Generating indices...")
combInd = Combinations(n, 4)
variations = Combinations(4).list()[1:] # 15 nontrivial combinations of 4 columns

A = matrix(G, 0, n)
# Generate A
print("Building matrix A. May take some time.")
for ind in combInd:
	for var in variations:
		newRow = matrix(G, 1, n)
		nInd = [ind[i] for i in var]
		newRow[:,nInd] = 1
		A = A.stack(newRow)

# row reduced echelon L, transformation T
print("Transforming.")
E = A.extended_echelon_form()
L, T = E[:, :n], E[:, n:]
Tinv = T.inverse()

print("Testing with gold function:")
xT = []		# top of columns
xB = []		# bottom of columns

for i in range(n):
	xT.append(kVecTok2field(vector(w^i), 0))
	xB.append(kVecTok2field(vector(f(w^i)),3))
xT = vector(G, xT)
xB = vector(G, xB)
x = xT + xB

inhom = A*x

print("Generating conditions")
VBasis = V.basis()
optionBasis = VBasis[:]
TB = T[n:, :]
for row in tqdm(TB.rows()):
	row = matrix(row)
	print("our option space had dimension %d." % len(optionBasis))
	zero = row.right_kernel()
	zeroBasis = zero.basis()
	optionBasis = [vec for vec in optionBasis if vec not in zeroBasis]
	print("now it has dim %d."% len(optionBasis))
	
#result = [0]

#while 0 in result:
	#bottom = W.random_element().list()
	#while 0 in bottom:
		#bottom = W.random_element().list()
	#right = []
	#for i in range(2^m-1):
		#col = vector(w^i).list() + vector(bottom[i]).list()
		#res = 0
		#for j in range(2*m):
			#print j
			#if col[j] != 0:
				#print("adding %s"%str(basis[j]))
				#res += basis[j]
		##print(col, res)
		#right.append(res)
	#right = vector(G, right)
	#result = A*right
	#print("just checked %s..."%str(right))
	
#print("%s did the job with lower %s."%(str(right), str(bottom)))	
