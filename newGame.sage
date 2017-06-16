import multiprocessing as mp
from time import time
import pickle
from numpy import random, array_split

m = 6
F = GF(2)
K.<w> = GF(2^m, 'w')
V = VectorSpace(F, 2*m) 		# linear independence checking

class Tree():
	def __init__(self, elem, parent=None, children=None):
		self.elem = elem
		self.parent = parent
		self.children = children
		
	def depth(self):
		depth = 0
		children = self.children
		while children != None:
			
			depth += 1
			candidates = []
			for child in children:
				if child.children != None:
					candidates += child.children
			if len(candidates) > 0:
				children = list(candidates)
			else:
				children = None
		return depth
		
def leafWorker(leaf, nWorkers = 8):
	used = [leaf.elem]
	parent = leaf.parent
	while parent != None:
		if parent.elem != K(0): used.append(parent.elem)
		parent = parent.parent
				
	# Generate Matrix up til now
	nCols = len(used)
	used.reverse()
	lower = matrix(GF(2), [vector(opt) for opt in used] ).transpose()
	upper = matrix(GF(2), [vector(w^i) for i in range(nCols)]).transpose()
	mat = upper.stack(lower)
						
	# Check possible options 
	options	= [elem for elem in K if elem not in used and elem != K(0)]
	print("Now investigating %d options to append to \n%s."%(len(options),mat.str()))
	children = []
	for option in options:
		if checkOption(mat, option, nWorkers):
			children.append(Tree(option, leaf))
	return children
			
def updateTreeMulti(tree, maxDepth = 3, nWorkers = mp.cpu_count(), leavesMax = floor(mp.cpu_count()/2)):
	depth = tree.depth()
	print("Tree depth is %d"%depth)
	while depth < maxDepth:
		if depth != 0:
			leaves = nodesOfRelDepth(tree, depth)
			
			if __name__ == "__main__":
				#p = mp.Pool(processes= nWorkers)
				if len(leaves) > leavesMax: leaves = random.choice(leaves, leavesMax).tolist()
				for leaf in leaves:
					leaf.children = leafWorker(leaf, nWorkers)
				#res = p.map(leafWorker, leaves)
				#for i in range(len(leaves)):
					#if len(res[i]) > 0:
						#leaves[i].children = res[i]
				
			newDepth = tree.depth()
			if newDepth != depth:
				depth = newDepth
			else:
				print("dead end")
				return False	
			
		else:
			# Intialization
			options = list(K.list())
			options.remove(K(0))
			tree.children = []
			for option in options:
				optionTree = Tree(option, tree)
				tree.children.append(optionTree)
			depth = tree.depth()
	return True
		
def updateTree(tree, maxDepth = 3):
	depth = tree.depth()
	print("Tree depth is %d"%depth)
	while depth < maxDepth:
		if depth != 0:
			leaves = nodesOfRelDepth(tree, depth)
			print("in first if")
			for leaf in leaves:
				used = [leaf.elem]
				parent = leaf.parent
				while parent != None:
					if parent.elem != K(0): used.append(parent.elem)
					parent = parent.parent
				
				# Generate Matrix up til now
				nCols = len(used)
				used.reverse()
				lower = matrix(GF(2), [vector(opt) for opt in used] ).transpose()
				upper = matrix(GF(2), [vector(w^i) for i in range(nCols)]).transpose()
				mat = upper.stack(lower)
								
				# Check possible options 
				options	= [elem for elem in K if elem not in used and elem != K(0)]
				print("Now investigating %d options to append to \n%s."%(len(options),mat.str()))

				children = []
				for option in options:
					if checkOption(mat, option):
						children.append(Tree(option, leaf))
				if len(children) > 0:
					leaf.children = children
			newDepth = tree.depth()
			if newDepth != depth:
				depth = newDepth
			else:
				print("dead end")
				return False
							
			
		else:
			# Intialization
			options = list(K.list())
			options.remove(K(0))
			tree.children = []
			for option in options:
				optionTree = Tree(option, tree)
				tree.children.append(optionTree)
			depth = tree.depth()
	return True

def nodesOfRelDepth(root, depth, currentDepth = 0):
	results = []
	if depth != 0:
		children = root.children 
		if children != None:
			if currentDepth + 1 == depth:
				results += children
			else:
				for child in children:
					results += nodesOfRelDepth(child, depth, currentDepth + 1)
		return results
	else:
		return [root]	
	
	
# Generate index lookup for checking linear independence of columns
# indices[nCols] = Combinations(nCols,4)
print("Generating indices for later...")
indices = dict()
for i in range(4, 2^m):
	indices[i] = Combinations(i, 4).list()
print("Done, proceeding.")

def checkOption(mat, opt, nWorkers): 		
	# vector to add to matrix
	N = mat.ncols()
	newColUpper = list(vector(w^N))
	newColLower = list(vector(opt))
	newCol = matrix(F, newColUpper+newColLower).transpose()
	
	newMat = mat.augment(newCol)
	res = checkRanksMulti(newMat, nWorkers)
	#if res:
		#print("\n%s\n works"%newMat.str())
	return res
	
def checkRanks(mat):
	N = mat.ncols()
	if N < 4:
		cols = mat.columns()
		if V.are_linearly_dependent(cols):
			return False
		else:
			return True
	else:
		inds = indices[N]
		for ind in inds:
			cols = mat[:, ind].columns()
			if V.are_linearly_dependent(cols):	# faster than rank checking it seems
				return False
		else:
			return True	

def checkRanksMulti(mat, nWorkers):
	N = mat.ncols()
	if N < 4:
		cols = mat.columns()
		if V.are_linearly_dependent(cols):
			return False
		else:
			return True
	else:
		inds = indices[N]
		chunks = array_split(inds, nWorkers)
		chunks = [chunk.tolist() for chunk in chunks]
		workers = []
		manager = mp.Manager()
		returnVals = manager.list()
		for chunk in chunks:
			p = mp.Process(target=checkIndices, args=(mat, chunk,returnVals))
			p.start()
			workers.append(p)
			
		while True:
			for worker in workers:
				if not worker.is_alive():
					workers.remove(worker)
					break
			if False in returnVals and len(workers) > 0:
				print("Terminating processes")
				for worker in workers:
					worker.terminate()
					workers.remove(worker)
				return False
			elif False not in returnVals and len(workers) == 0:
				return True

def checkIndices(mat, indexList, retVals):
	for ind in indexList:
		cols = mat[:, ind].columns()
		if V.are_linearly_dependent(cols):	# faster than rank checking it seems
			retVals.append(False)
			return 
	else:
		return 

	
root = Tree(K(0))
nCols = 0
maxDepth = 3
newRoot = root

def chooseNewRoot(parent, maxDepth):
	if parent.children != []:
		nMaxBranches = []
		for child in parent.children:
			nMaxBranches.append(len(nodesOfRelDepth(child, maxDepth - 1)))
		newRoot = parent.children[ nMaxBranches.index(max(nMaxBranches)) ]
		
		return newRoot
	else:
		return False
	
while nCols < 2^m - 1:
	if updateTreeMulti(newRoot, maxDepth, 9):
		## Select new root among children
		newRoot = chooseNewRoot(newRoot, maxDepth)
		
		if newRoot:
			if nCols == 0:
				nCols = maxDepth
			else:
				nCols += 1
		else:
			print("No options left on first stage.")
			break
	else:
		while newRoot.parent != None:
			if newRoot.parent.children != None:
				newRoot.parent.children.remove(newRoot)
				newestRoot = chooseNewRoot(newRoot.parent, maxDepth)
				if newestRoot:
					newRoot.parent.children.append(newRoot)
					newRoot = newestRoot
					break
				else:
					newRoot = newRoot.parent
					nCols -= 1
			else:
				newRoot = newRoot.parent
				nCols -= 1
			
## print a Solution
sols = nodesOfRelDepth(root, 2^m - 1)
print(len(sols))
if len(sols) > 0:
	sol = sols[0]

	used = [sol.elem]
	parent = sol.parent
	while parent != None:
		if parent.elem != K(0): used.append(parent.elem)
		parent = parent.parent
					
	## Generate Matrix up til now
	nCols = len(used)
	used.reverse()
	lower = matrix(GF(2), [vector(opt) for opt in used] ).transpose()
	upper = matrix(GF(2), [vector(w^i) for i in range(nCols)]).transpose()
	mat = upper.stack(lower)
									
	print("A solution is \n%s."%(mat.str()))
	with open("TreeAPN%.2f"%time(), "w") as f:
		f.write(mat.str())
with open("Tree","w") as f:
	pickle.dump(root, f)
	
