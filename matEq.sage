from sage.geometry.hyperplane_arrangement.affine_subspace import AffineSubspace
from tqdm import tqdm
import pickle
import os

cont = False

for k in range(1,6):
	print("k is %d"%k)
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

	if not os.path.isfile("MatrixADim%d.data"%k):
		A = matrix(G, 0, n)
		# Generate A
		print("Building matrix A. May take some time.")
		for ind in tqdm(combInd):
			for var in variations:
				newRow = matrix(G, 1, n)
				nInd = [ind[i] for i in var]
				newRow[:,nInd] = 1
				A = A.stack(newRow)
		with open("MatrixADim%d.data"%k,"w") as f:
			pickle.dump(A, f)
	else:
		with open("MatrixADim%d.data"%k,"w") as f:
			A = pickle.load(f) 

	# row reduced echelon L, transformation T
	print("Transforming.")
	if not os.path.isfile("LTDim%d.data"%k):
		E = A.extended_echelon_form()
		L, T = E[:, :n], E[:, n:]
		Tinv = T.inverse()
		with open("LTDim%d.data"%k, "w") as f:
			pickle.dump( (L,T,Tinv) , f)
	else:
		with open("LTDim%d.data"%k, "r") as f: 
			(L,T,Tinv) = pickle.load(f)
	
if cont:
	print("Testing with gold function:")
	xT = []		# top of columns
	xB = []		# bottom of columns

	for i in range(n):
		xT.append(kVecTok2field(vector(w^i), 0))
		xB.append(kVecTok2field(vector(f(w^i)),k))
	xT = vector(G, xT)
	xB = vector(G, xB)
	x = xT + xB

	inhom = A*x

	print("Generating conditions")
	TRed = T[:, :n]
	TinvRed = Tinv[:, :n]
	xi = TinvRed*xT
	avoid = set()
	for i in range(m):
		row = TinvRed[i, :]
		comp = matrix(G, 1,1, xi[i])
		rowExt = row.augment(comp)
		print comp
		if row.rank() == rowExt.rank():
			kern = row.right_kernel()
			sol = row.solve_right(comp).columns()[0]
			print sol
			solSpace = AffineSubspace(sol,kern)
			avoid.add(solSpace)
			print solSpace
			print("sol space %d."%i)
		else:
			print("Not solvable")

	print("Testing bottom options")
	sub = []
	for i in range(n):
		sub.append(kVecTok2field(vector(w^i), k))
		
	variate = [xB]
	solutions = []
	while len(solutions) < 1:
		for var in variate:
			count = 0
			print("Variate is %s."%str(var))
			while count < n-1:
				cp = var[:]
				for elem in sub:
					cp[count] = elem
					print("Testing \n%s"%str(cp))
					for affSpace in avoid:
						if cp in affSpace:
							print("no luck %d "%len(solutions))
							break
					else:
						print("sol!!!!")
						solutions.append(cp[:])
					if not cp in variate: variate.append(cp)	
				count += 1
		

	
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
