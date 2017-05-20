from multiprocessing import Process, Manager 
from time import time
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
			
	if nCols < 2^m-2:
		for col, opt in lookUp:
			Mat = mat.augment(col)
			print("While having %d sols, investigating \n%s\n"%(len(solList),Mat.str()))
			nextCol(Mat, opt, solList)
	else:
		for col,opt in lookUp:
			Mat =mat.augment(col)
			print("""SUCCESS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			------------------------------------------------------------------------------------------------------------
			%s
			------------------------------------------------------------------------------------------------------------"""%((Mat.str())))
			solList.append(Mat)
			with open("sols_%.2f.data"%time(), "w+") as f:
				pickle.dump(solList,f)
def search(workers = 2):
	if workers != None:
		
		unused = V[1:]
		workerObj = list()
		if __name__ =='__main__':
			manager = Manager()
			sols = manager.list()
			for i in range(0,workers):
				elem = unused[i]
				lowerCol = matrix( [elem] ).transpose()
				mat = upper[:,0].stack(lowerCol)
				unusedC = list(unused)
				unusedC.remove(elem)
				p = Process(target=nextCol, args=(mat, unusedC, sols))
				p.start()
				workerObj.append(p)
		i = workers 
		while len(workerObj)>0:
			for p in workerObj:
				if p.exitcode != None and len(unused)>0:
					i+=1
					workerObj.remove(p)
					elem = unused[i]
					unusedC = list(unused)
					lowerCol = matrix( [elem] ).transpose()
					mat = upper[:,0].stack(lowerCol)
					unusedC.remove(elem)
					p = Process(target=nextCol, args=(mat, unusedC, sols))
					p.start()
					workerObj.append(p)
				elif p.exitcode != None:
					workerObj.remove(p)
					
		return 	list(sols)
		
	else:
		unused = V[1:]
		sols = list()
		for elem in unused:
			lowerCol = matrix([elem]).transpose()
			mat = upper[:,0].stack(lowerCol)
			unusedC = list(unused)
			unusedC.remove(elem)
			nextCol(mat, unusedC, sols)
			
Sols = search(None)

with open('sols.data','w+') as f:
	pickle.dump(Sols,f)
