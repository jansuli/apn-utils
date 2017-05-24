class APN():
	def __init__(self, domain, codomain, expression, primName = 'a'):
		self.domain = domain
		self.codomain = codomain
		self.expression = expression
		
		self.domainVecSpace = domain.vector_space()
		self.codomainVecSpace = codomain.vector_space()
		
		self.primName = primName
		self.primElem = domain.primitive_element()
	def func(self,x):	
		return sage_eval(self.expression,locals = {'x':x, self.primName: self.primElem} )
		
	def vecFunc(self, vec):
		fieldElem = self.domain.zero()
		i = 0
		for coordinate in vec:
			fieldElem += coordinate*self.primElem^i 
			i += 1
		return vector(self.func(fieldElem))
		
	def componentFunctions(self):
		components = dict()
		def component(b):
			def Fb(v):
				res = b.dot_product(self.vecFunc(v))
				return res
			return Fb
			
		for b in self.codomainVecSpace:
			b.set_immutable()
			components[b] = component(b)
		self.components = components
