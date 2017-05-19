try:
	tqdm = tqdm
	Pool = Pool
	partial = partial
	time = time
	np = numpy
except NameError:
	import numpy as np
	from tqdm import tqdm
	from multiprocessing import Pool
	from functools import partial
	import time
	
	
def checkIndependentCols(mat,D,indexList):
	r = mat[:, indexList].rank()
	if r != D:
		return False 
				
def searchDependentCols(mat,D,indexList):
	r = mat[:, indexList].rank()
	if r == D:
		return True

def minimumDistance(mat,d, workers=2, chunks=4):
	''' Taking a the parity check matrix 'mat' of a Code this function returns True if the Code has minimum distance d, else False.
	
	It does so by checking each tuple of 4 different columns for linear independce (via matrix rank) and looking for a 5-tuple of linear dependent columns.
	This process can get quite tedious as for the first step alone we have to check binomial{mat.ncols()}{4} combinations.'''
	
	n = mat.ncols()
	
	if workers == 1:
		indices = Combinations(range(n), d-1).list()
		D = d-1
		print("Looking for indipendent %d columns."%(D))
		for ind in tqdm(indices):
			r = mat[:, ind].rank()
			
			if r != D:
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

	else:	
		try:
			if __name__ == '__main__':
				indices = Combinations(range(n), d-1).list()
				indices = np.array_split(indices, workers*chunks)
				indices = [arr.tolist() for arr in indices]
				for indexGroup in tqdm(indices):
					pool = Pool(processes = workers)
					func = partial(checkIndependentCols, mat, d-1)
					out = pool.map(func, indexGroup)
					pool.close()
					if False in out:
						return False
				else:
					indices = Combinations(range(n), d).list()
					indices = np.array_split(indices, workers*chunks)
					indices = [arr.tolist() for arr in indices]
					for indexGroup in tqdm(indices):
						pool = Pool(processes = workers)
						func = partial(searchDependentCols, mat, d)
						out = pool.map(func, indexGroup)
						pool.close()
						if True in out:
							return True
					else:
						return False
		except KeyboardInterrupt:
			pool.terminate()
			pool.close()
		
	
