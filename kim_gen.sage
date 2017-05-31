from tqdm import tqdm
from time import sleep
from multiprocessing import Pool, Manager
from numpy import array_split
from functools import partial
import pickle
m = 10
q = 2^(m/2)

K.<w> = GF(2^m, 'w')
L = [K(0)] + [e for e in K if e!= 0 and (q-1) % e.multiplicative_order() == 0] # subfield
k = list(K)
k.remove(K(0))

print L

cands = []
for beta in tqdm(k):
	for gamma in k:
		for scalar in L:
			if beta == scalar*gamma:
				break
		else:
			cands.append( (gamma,beta) )

kappa = 1

def checkPair(sols, pair):
	cont= True
	gamma = pair[0]
	beta = pair[1]
	for a in K:
		if cont:
			#print(a)
			for b in K:
				if cont:
					for c in K:
						# Checking x maps to tr^m_k(beta, k(x,y) bijectivity
						# precheck
						pre = beta + a*beta + b*beta + c*beta
						#tqdm.write(str(pre))
						if not pre in L:
							left = (beta+a*beta+b*beta+c*beta+beta^q+a^q*beta^q+b^q*beta^q+c^q*beta^q)*(w*beta + a*w^q*beta + b*w*beta +c*w^q*beta + w^q*beta^q+a^q*w*beta^q + b^q*w^q*beta^q + c^q*w*beta^q)^(2^kappa)
							right = (beta + a*beta + b*beta + c*beta + beta^q + a^q*beta^q + b^q*beta^q + c^q*beta^q)^(2^kappa) * (w^(2^kappa)*beta +a*w^(2^kappa)*beta+b*w^(2^kappa*q)*beta +c*w^(2^kappa *q)*beta + w^(2^kappa *q)*beta^q + a^q *w^(2^kappa *q)*beta^q + b^q *w^(2^kappa)*beta^q + c^q *w^(2^kappa)*beta^q)
							
							if left == right:
								# check other way round
								# precheck
								pre=(w^(2^kappa + 1))*gamma + a*(w^(2^kappa + q))*gamma + b*(w^(2^kappa*q+1))*gamma + c*(w^(2^kappa*q+q))*gamma
								if not pre in L:
									#tqdm.write("gamma precheck success")
									left= (w^(2^kappa +1)*gamma + a*w^(2^kappa + q)*gamma + b*w^(2^kappa*q+1)*gamma + c*w^(2^kappa*q+q)*gamma + w^(2^kappa*q+q)*gamma^q + a^q*w^(2^kappa*q+1)*gamma^q + b^q*w^(2^kappa +q) *gamma^q + c^q * w^(2^kappa+1) *gamma^q)*(w^(2^kappa)*gamma + a*w^(2^kappa)*gamma + b*w^(2^kappa*q)*gamma + c*w^(2^kappa*q)*gamma + w^(2^kappa*q)*gamma^q +a^q *w^(2^kappa*q)*gamma^q + b^q*w^(2^kappa)*gamma^q + c^q*w^(2^kappa)*gamma^q)^(2^kappa)
									right = (w^(2^kappa + 1)*gamma + a*w^(2^kappa + q) *gamma + b*w^(2^kappa * q +1)*gamma +c*w^(2^kappa *q+q)*gamma +w^(2^kappa*q+q)*gamma^q + a^q *w^(2^kappa*q+1)*gamma^q + b^q*w^(2^kappa+q)*gamma^q + c^q*w^(2^kappa +1)*gamma^q)^(2^kappa) * ( w*gamma + a*w^q + b*w*gamma + c*w^q *gamma+ w^q*gamma^q + a^q*w^q*gamma^q + b^q*w^q*gamma^q + c^q*w*gamma^q ) 
									
									if left == right:
										tqdm.write("succes ----------------------------------")
										tqdm.write( str((gamma, beta, a,b,c)))
										sols.append( (gamma,beta,a,b,c) )
										print(len(sols))
										cont = False
										break
i = 0			
sols = []
while len(sols)<1:
	checkPair(sols,cands[i])
	i+=1

sol=sols[0]

genericComponent = "(delta + a*delta + b*delta + c * delta + delta^q + a^q*delta^q + b^q*delta^q +c^q*delta^q)*x^(2^kappa+1) +(w*delta +a*w^q*delta+b*w*delta+c*w^q*delta+w^q*delta^q + a^q*w*delta^q + b^q*w^q*delta^q +c^q*w*delta^q)*x^(2^kappa)*y + (w^(2^kappa)*delta+a*w^(2^kappa)*delta+b*w^(2^kappa)*delta+c*w^(2^kappa)*delta+ w^(2^kappa*q)*delta^q+a^q*w^(2^kappa*q)*delta^q+b^q*w^(2^kappa)*delta^q + c^q*w^(2^kappa)*delta^q)*x*y^(2^kappa)+(w^(2^kappa +1)*delta +a*w^(2^kappa + q)*delta + b*w^(2^kappa*q+1)*delta +c*w^(2^kappa*q+q)*delta + w^(2^kappa*q+q)*delta^q + a^q*w^(2^kappa*q+1)*delta^q + b^q*w^(2^kappa+q)*delta^q + c^q *w^(2^kappa+1)*delta^q)*y^(2^kappa +1)"

def returnFunction(sol, pos):
	'''0 -> gamma, 1 -> beta'''
	def func((x,y)):
		res = sage_eval(genericComponent, locals={'delta':sol[pos],'a':sol[2],'b':sol[3], 'c':sol[4], 'x':x, 'y':y, 'w':w, 'q':q, 'kappa':kappa})
		return res
	return func
	
F = returnFunction(sol,0)
G = returnFunction(sol,1)
def returnComp(func,pos):
	def F((x,y)):
		#print ("Handling %s with return comp: %s."%(str((x,y)),(x,y)[pos])) 
		#print func((x,y))
		#print ("\n")
		return ((x,y)[pos], func((x,y)))
	return F

F1 = returnComp(F,0)
F2 = returnComp(G,1)

F1Table = dict()
F2Table = dict()
LtimesL = []
for a in L:
	for b in L:
		LtimesL.append( (a,b) )

for vec in tqdm(LtimesL):
	res = F1(vec)
	if res not in F1Table.values():F1Table[vec] = res		
for vec in tqdm(LtimesL):
	res = F2(vec)
	if res not in F2Table.values():F2Table[vec] = res
		
print len(F1Table)
print len(F2Table)

F1TableInv = dict()
for k,v in F1Table.iteritems():
	F1TableInv[v] = k

	
for vec in LtimesL:
	print("%s maps to:"%str(vec))
	intermediate = F1TableInv[vec]
	res = F2Table[intermediate]
	print("%s \n"%str(res))
