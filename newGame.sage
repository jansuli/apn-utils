from tqdm import tqdm

m = 6
F = GF(2)
K.<w> = GF(2^m, 'w')
V = VectorSpace(F, 2*m) 		# linear independence checking

class Tree():
	def __init__(self, elem, parent=None, children=None):
		self.elem = elem
		
		if children != None:
			self.children = children
		else:
			self.children = []
			
		if not parent:
			self.parent = None
		else:
			self.parent = parent

# Generate index lookup for checking linear independence of columns
# indices[nCols] = Combinations(nCols,4)
print("Generating indices for later...")
indices = dict()
for i in range(4, 2^m):
	indices[i] = Combinations(i, 4)
print("Done, proceeding.")

def checkOption(mat, opt): 		
	# vector to add to matrix
	N = mat.ncols()
	newColUpper = list(vector(w^N))
	newColLower = list(vector(opt))
	newCol = matrix(F, newColUpper+newColLower).transpose()
	
	newMat = mat.augment(newCol)
	
	return checkRanks(newMat)
	
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
			
def chooseNext(trees):
	lengths = [len(tree.children) for tree in trees]
	most = max(lengths)
	
	if len(trees) != lengths.count(most):
		return trees[ lengths.index(most) ]
	else:
		return trees[randint(0,len(trees)-1)] 
			
# Initialization
trees = list()
for i in tqdm(range(0,2^m-1), desc="Initiator"):
	elem = w^i
	
	possibleOptions = copy(K.list())
	if K(0) in possibleOptions:
		possibleOptions.remove(K(0))
	
	for option in tqdm(possibleOptions):
		elemTree = Tree(option)
		children = []
		optionsLeft = copy(possibleOptions)
		optionsLeft.remove(elemTree.elem)
		
		initialUpper = vector(elem).list()
		initialLower = vector(option).list()
		initialMatrix = matrix(F, initialUpper+initialLower).transpose()
		
		for opt in optionsLeft:
			if checkOption(initialMatrix, opt):
				optTree = Tree(opt,parent=elemTree)
				children.append(optTree)
	
		elemTree.children = children
		trees.append(elemTree)

N = 0
while N < 2^m-1:
	tree = chooseNext(trees)
	trees = []
	
	# construct matrix until now
	N = 1
	lower = [vector(tree.elem)]
	parent = tree.parent
	while parent != None:
		N += 1
		lower.append(vector(parent.elem))
		parent = parent.parent
	lower.reverse()
	upper = [vector(w^i) for i in range(N)]
	
	upperMat = matrix(F, upper).transpose()
	lowerMat = matrix(F, lower).transpose()
	initialMatrix = upperMat.stack(lowerMat)
	
	print(initialMatrix.str())
	
	for optionTree in tqdm(tree.children, desc="Options"):
		option = optionTree.elem
		newCol = matrix(F, vector(w^N).list() + vector(option).list()).transpose()
		testMatrix = initialMatrix.augment(newCol)
		optionsLeft = [elem for elem in K if elem != 0 and elem not in lower + [option]]
		for opt in optionsLeft:
			if checkOption(testMatrix,opt):
				newTree=Tree(opt,parent=optionTree)
				optionTree.children.append(newTree)
		trees.append(optionTree)
			
