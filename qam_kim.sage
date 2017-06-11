from tqdm import tqdm
import pickle
from time import time

m=10
q=2^(m/2)

K.<w> = GF(2^m, 'w')
L = [K(0)] + [e for e in K if e!= 0 and (q-1) % e.multiplicative_order() == 0] # subfield
LtL = CartesianProduct(L,L)

LtL2K = dict()
K2LtL = dict()
### Isomorphism:
for vec in LtL:
	elem = vec[0] + w*vec[1]
	LtL2K[vec] = elem
for key, val in LtL2K.iteritems():
	K2LtL[val] = key

### Functions
def uniKim( (a,b,c) ):
	def kim(x):
		return sage_eval("x^3 + a*x^32 + b*x^65 + c*x^96", locals={'a':a,'b':b,'c':c, 'x':x})
	return kim
	
def bivKim( (a,b,c) ):
	ev = "(1+a+b+c)*x^(2^kappa +1) + (w+a*w^q +b*w +c*w^q)*x^(2^kappa)*y + (w^(2^kappa) + a*w^(2^kappa) + b*w^(2^kappa*q) + c*w^(2^kappa *q))*x*y^(2^kappa) + (w^(2^kappa + 1) +a*w^(2^kappa + q) +b*w^(2^kappa*q +1) +c*w^(2^kappa*q+q))*y^(2^kappa + 1)"
	#ev = "(1+a+b+c)*x^3 + (w+a*w^32+b*w+c*w^32)*x^2*y + (w^2+a*w^2 +b^64 +c*w^64)*x*y^2 + (w^3 + a*w^34 +b*w^65 + c*w^96)*y^3"
	def kim((x,y)):
		return sage_eval(ev, locals={'a':a,'b':b,'c':c, 'w':w, 'x':x,'y':y, 'q':q, 'kappa':1})
	return kim
	
def traceKim( (a,b,c), delta):
	ev = "(delta + a*delta + b*delta + c * delta + delta^q + a^q*delta^q + b^q*delta^q +c^q*delta^q)*x^(2^kappa+1) +(w*delta +a*w^q*delta+b*w*delta+c*w^q*delta+w^q*delta^q + a^q*w*delta^q + b^q*w^q*delta^q +c^q*w*delta^q)*x^(2^kappa)*y + (w^(2^kappa)*delta+a*w^(2^kappa)*delta+b*w^(2^kappa*q)*delta+c*w^(2^kappa*q)*delta+ w^(2^kappa*q)*delta^q+a^q*w^(2^kappa*q)*delta^q+b^q*w^(2^kappa)*delta^q + c^q*w^(2^kappa)*delta^q)*x*y^(2^kappa)+(w^(2^kappa +1)*delta +a*w^(2^kappa + q)*delta + b*w^(2^kappa*q+1)*delta +c*w^(2^kappa*q+q)*delta + w^(2^kappa*q+q)*delta^q + a^q*w^(2^kappa*q+1)*delta^q + b^q*w^(2^kappa+q)*delta^q + c^q *w^(2^kappa+1)*delta^q)*y^(2^kappa +1)"
	#ev = "(delta + a*delta + b*delta + c * delta + delta^32 + a^32*delta^32 + b^32*delta^32 +c^32*delta^32)*x^3 +(w*delta +a*w^32*delta+b*w*delta+c*w^32*delta+w^32*delta^32 + a^32*w*delta^32 + b^32*w^32*delta^32 +c^32*w*delta^32)*x^2*y + (w^2*delta+a*w^2*delta+b*w^(64)*delta+c*w^(64)*delta+ w^64*delta^32+a^32*w^(64)*delta^32+b^32*w^(2)*delta^32 + c^32*w^(2)*delta^32)*x*y^2+(w^3*delta +a*w^34*delta + b*w^65*delta +c*w^96*delta + w^96*delta^32 + a^32*w^65*delta^32 + b^32*w^34*delta^32 + c^32 *w^3*delta^32)*y^3"
	def tr((x,y)):
		return sage_eval(ev, locals={'a':a,'b':b,'c':c, 'w':w, 'x':x, 'y':y, 'delta':delta, 'q':q, 'kappa':1})
	return tr
	
### Find Gamma and Beta:
def findGammaBeta():
	a = K.random_element()
	b = K.random_element()
	c = K.random_element()

	bKim = bivKim((a,b,c))

	for gamma in tqdm(K):
		F = traceKim((a,b,c),gamma)
		for vec in LtL:
			res = F(vec)
			if res != K2LtL[bKim(vec)][0]:
				break
		else:
			print("possible gamma is %s"%str(gamma))
			for beta in K:
				G = traceKim((a,b,c), beta)
				for vec in LtL:
					if G(vec) != K2LtL[bKim(vec)][1]:
						break
				else:
					print("Gamma and Beta are:")
					print(gamma,beta)
					return gamma,beta
(gamma, beta) = (w^8 + w^6 + w^5 + w^4 + w^2 + 1, w^7 + w^5 + w^4 + w^3 + w)

### CCZ permutation condition checks:
def checkParams( (a,b,c), gamma, beta):
	kappa = 1
	factor1 = (beta+a*beta+b*beta+c*beta+beta^q+a^q*beta^q+b^q*beta^q+c^q*beta^q)
	if factor1 != 0:
		#tqdm.write("first check successfull")
		left = factor1*(w*beta + a*w^q*beta + b*w*beta +c*w^q*beta + w^q*beta^q+a^q*w*beta^q + b^q*w^q*beta^q + c^q*w*beta^q)^2
		right = factor1^2 * (w^2*beta +a*w^2*beta+b*w^64*beta +c*w^64*beta + w^64*beta^q + a^q *w^64*beta^q + b^q *w^2*beta^q + c^q *w^2*beta^q)
		
		if left == right:
			#tqdm.write("second check successfull")
			#pre=(w^(2^kappa + 1))*gamma + a*(w^(2^kappa + q))*gamma + b*(w^(2^kappa*q+1))*gamma + c*(w^(2^kappa*q+q))*gamma
			factor2 = (w^3*gamma + a*w^34*gamma + b*w^65*gamma + c*w^96*gamma + w^96*gamma^q + a^q*w^65*gamma^q + b^q*w^34 *gamma^q + c^q * w^3 *gamma^q)
			if factor2 != 0:
				#tqdm.write("third check successfull")
				left= factor2*(w^2*gamma + a*w^2*gamma + b*w^64*gamma + c*w^64*gamma + w^64*gamma^q +a^q *w^64*gamma^q + b^q*w^2*gamma^q + c^q*w^2*gamma^q)^2
				right = factor2^2 * ( w*gamma + a*w^q + b*w*gamma + c*w^q *gamma+ w^q*gamma^q + a^q*w^q*gamma^q + b^q*w^q*gamma^q + c^q*w*gamma^q ) 
				
				if left == right:
					return True
	return False


#### QAM - Checking ###
M = matrix(K, m)
for i in range(10):
	for j in range(10):
		M[i,j] = (w^j)^(2^i)
		
indices = Combinations(m).list()[1:]

def checkQAM( (a,b,c) ):
	
	C = matrix(K, m)

	C[5,1] = C[1,5] = a
	C[6,0] = C[0,6] = b
	C[6,5] = C[5,6] = c

	H = M.transpose() * C * M
	
	for ind in indices:
		cols = H[:, ind].columns()
		linComb = vector(sum(cols))
		rCheck = matrix(vector(linComb[0]))
		for elem in linComb[1:]:
			rCheck = rCheck.stack(vector(elem))
		
		if rCheck.rank() != m-1:
			return False
	else:
		return True
		
#a = K.random_element()
#print("a ist %s"%str(a))
apns = list()
try:
	for a in tqdm(K, desc="a"):
		for b in tqdm(K, decs="b"):
			for c in tqdm(K,desc='c'):
				params = checkParams((a,b,c),gamma,beta) 
				if params:
					tqdm.write("%s fit the conditions..."%(str((a,b,c))))
					apns.append((a,b,c))
					if checkQAM((a,b,c)):
						tqdm.write("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Equiv to a permutation !!!!!!!!!!!!!!!!!!!!!!!!!")
						with open("PERMUTATION%f.DATA"%time(),"w") as f:
							pickle.dump( (a,b,c), f)
					else:
						tqdm.write("sadly not apn...")
except KeyboardInterrupt:
	with open("apns","w") as f:
		pickle.dump(apns, f)


#cond = False
#while not cond:
	#a = K.random_element()
	#b = K.random_element()
	#c = K.random_element()
	#print( a,b,c )
	#preCond = checkQAM( (a,b,c) )
	#if preCond:
		#cond = checkParams((a,b,c),gamma,beta)
	#print preCond,cond
