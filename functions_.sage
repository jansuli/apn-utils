def alternatingBilinearForm(func, vecSpace):
	'''Maps boolean function 'func':F^m -> F (we will later use the component functions) to alternating bilinear form.'''
	
	def Bf(x,y):
		return func(x+y) + func(x) + func(y) + func(vecSpace.zero())
	return Bf
	
def alternatingMatrix(form, basis):
	size = basis[0].parent().dimension()
	mat = matrix(GF(2), size)
	for i in range(1,size):
		for j in range(0,i):
			mat[i,j] = form(basis[i], basis[j])
	mat = mat + mat.transpose()
	return mat
