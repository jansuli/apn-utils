from multiprocessing import Process, Manager 
import pickle
m=6
K.<a> = GF(2^m, 'a')
V = K.vector_space()
BigV = VectorSpace(GF(2), m*2)
print("The minimal polynomial of 'a' is:")
print(a.minimal_polynomial())
k = K.list()
upper = matrix(k[1:])


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

upper = GFtoBinMatrix(upper,m)
manager = Manager()
sols = manager.list()

def nextCol(mat, options,solList):
	nCols = mat.ncols()
	if nCols < 3:
		combIndices = Combinations(range(nCols+1))
	else:
		combIndices = Combinations(range(nCols+1), 4)
	combIndices = [ind for ind in combIndices if nCols in ind and len(ind) > 1]
	# generate Lookup with apn-requirement
	columnSpaceBase = list()
	for ind in combIndices:
		cols = mat[:,ind[:-1]].columns()
		columnSpaceBase += cols
	columnSpace = BigV.subspace(columnSpaceBase) 
	lookUp = list()
	upperNew = upper[:,nCols]
	for option in options:
		lowerNew = matrix( [option]).transpose()
		colNew = upperNew.stack(lowerNew)
		if colNew not in columnSpace:
			optionsCopy = list(options)
			optionsCopy.remove(option)
			lookUp.append( (colNew, optionsCopy) )
		else:
			print ("%s was in linear dependent.")
			
	if nCols < 2^m - 2:
		for col, opt in lookUp:
			Mat = mat.augment(col)
			print("While having %d sols investigating \n%s\n"%(len(solList),Mat.str()))
			nextCol(Mat, opt, solList)
	else:
		for col,opt in lookUp:
			print("SUCCESS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
			sols.append(mat.augment(col))
def search(workers = 2):
	unused = V[1:]
	workerObj = list()
	if __name__ =='__main__':
		for i in range(0,workers):
			elem = unused[i]
			lowerCol = matrix( [elem] ).transpose()
			mat = upper[:,0].stack(lowerCol)
			unused.remove(elem)
			p = Process(target=nextCol, args=(mat, unused, sols))
			p.start()
			workerObj.append(p)

	while len(workerObj)>0:
		for p in workerObj:
			if p.exitcode != None and len(unused)>0:
				workerObj.remove(p)
				p = Process(target=nextCol, args=(mat,unused,sols))
				elem = unused[0]
				lowerCol = matrix( [elem] ).transpose()
				mat = upper[:,0].stack(lowerCol)
				unused.remove(elem)
				p = Process(target=nextCol, args=(mat, unused, sols))
				p.start()
				workerObj.append(p)
			elif p.exitcode != None:
				workerObj.remove(p)
				
	return 	list(sols)

with open('sols.data','w+') as f:
	pickle.dump(sols,f)
