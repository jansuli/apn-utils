from funcs import functions		# List of functions
from sage.coding.codecan.autgroup_can_label import LinearCodeAutGroupCanLabel

functionStrings = functions.replace("\n","").split(',')
functionStrings.remove('x^57')

n = 8
K.<g> = GF(2^n, 'g', repr="log")
KList = K.list()[1:]

def myFx(x):
	return g^250*x^192 + g^236*x^160 + g^182*x^144 + g^82*x^136 + g^82*x^132 + g^152*x^130 + g^132*x^129 + g^244*x^96 + g^201*x^80 + g^152*x^72 + g^131*x^68 + g^191*x^66 + g^63*x^65 + g^84*x^48 + g^156*x^40 + g^109*x^36 + g^221*x^34 + g^93*x^33 + g^239*x^24 + g^148*x^20 + g^97*x^18 + g^162*x^17 + g^91*x^12 + g^176*x^10 + g^254*x^9 + g^183*x^6 + g^188*x^5 + g^68*x^3

def builtGenerator(func):
	A = matrix(GF(2), 2*n, 0)
	for x in KList:
		top = matrix(vector(x)).transpose()
		bottom = matrix(vector(func(x))).transpose()
		col = top.stack(bottom)
		A = A.augment(col)
	return A

myMat = builtGenerator(myFx)
print ("Built matrix \n%s..."%myMat[:,:10].str())
C = LinearCode(myMat)

CAN = LinearCodeAutGroupCanLabel(C, algorithm_type="linear").get_canonical_form()

for fx in functionStrings:
	def f(x):
		return sage_eval(fx, locals = {'x':x, 'g':g})
		
	A = builtGenerator(f)
	print("Testing \n%s..."%A[:,:10].str())
	C2 = LinearCode(a)
	can = LinearCodeAutGroupCanLabel(C2, algorithm_type="linear").get_canonical_form()
	
	if CAN == can:
		print("APN CCZ-equivalent to\n%s\n"%fx)
		break
else:
	print("Code inequivalent to all tested functions :)")
	
	
