from functools import partial
import multiprocessing as mp
from time import time

m = 6
K.<w> = GF(2^m, 'w')
V = VectorSpace(GF(2), 2*m)

def f(x):
	return x^3 + x^10 + w*x^24

print("Analysing f(x).")
image = dict()
for elem in K:
	if elem != K(0):
		res = f(elem)
		if res in image:
			image[res].append(elem)
		else:
			image[res] = [elem]
print("f maps K* to %d different elements."%len(image))

print("Pre-Building our matrix")
k = K.list()[1:]
matK = matrix(K, 2, 0)
for upper in k:
	col = matrix(K,[upper, f(upper)]).transpose()
	matK = matK.augment(col)
print("Done, now generating indices.")

indices = Combinations(2^m - 1, 4).list()
indices = [ind for ind in indices if 2^m-2 in ind]

def checkOption(mat, upper, opt):
	colUpper = matrix(GF(2), vector(upper)).transpose()
	colLower = matrix(GF(2), vector(opt)).transpose()
	col = colUpper.stack(colLower)
	testMat = mat.augment(col)
	
	for ind in indices:
		cols = testMat[:, ind].columns()
		if V.are_linearly_dependent(cols):
			return
	else:
		return opt
		
checked = []
while len(checked) < 2^m-1:
	index = choice(range(2^m-1))
	top = mat[0, index]
	if not top in checked:
		bottom = matK[1, index]
		options = [elem for elem in K if elem != K(0) and elem != bottom]
		
		matKrem = matK[:, :index].augment(matK[:, index+1:])
				
		matGF = matrix(GF(2), 2*m, 0)
		for col in matKrem.columns():
			colGFUp = matrix(GF(2), vector(col[0])).transpose()
			colGFLo = matrix(GF(2), vector(col[1])).transpose()
			colGF = colGFUp.stack(colGFLo)
			matGF = matGF.augment(colGF)
		
		optionChecker = partial(checkOption, matGF, top)
		
		print("Investigating options for removed column %s.\n"%str(matK[:, index]))
		
		if __name__ == "__main__":
			p = mp.Pool(processes = mp.cpu_count())
			res = p.map(optionChecker, options)
			res = [r for r in res if r != None]
			p.close()
			p.join()
					
			if len(res)>0:
				for r in res:
					print("We found a new APN, Saving it")
					colGFUp = matrix(GF(2), vector(upper)).transpose()
					colGFLo = matrix(GF(2), vector(r)).transpose()
					col = colGFUp.stack(colGFLo)
						
					save = matGF.augment(col)
					with open("results/apn_%.1f"%time(),"w") as f:
						f.write(save.str())
						
		checked.append(top)
			
			
		

