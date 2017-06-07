from tqdm import tqdm
from time import sleep
from multiprocessing import Pool, Manager
from numpy import array_split
from functools import partial
import pickle
m = 10
q = 2^(m/2)
kappa = 1
K.<w> = GF(2^m, 'w')
L = [K(0)] + [e for e in K if e!= 0 and (q-1) % e.multiplicative_order() == 0] # subfield
k = list(K)
k.remove(K(0))		# for nontrivial trace functions

LtimesL = []
for a in L:
	for b in L:
		LtimesL.append( (a,b) )
		
LtL2K = dict()
for vec in LtimesL:
	LtL2K[vec] = vec[0] + w*vec[1]
K2LtL = dict()
for key,value in LtL2K.iteritems():
	K2LtL[value] = key

genericComponent = "(delta + a*delta + b*delta + c * delta + delta^q + a^q*delta^q + b^q*delta^q +c^q*delta^q)*x^(2^kappa+1) +(w*delta +a*w^q*delta+b*w*delta+c*w^q*delta+w^q*delta^q + a^q*w*delta^q + b^q*w^q*delta^q +c^q*w*delta^q)*x^(2^kappa)*y + (w^(2^kappa)*delta+a*w^(2^kappa)*delta+b*w^(2^kappa*q)*delta+c*w^(2^kappa*q)*delta+ w^(2^kappa*q)*delta^q+a^q*w^(2^kappa*q)*delta^q+b^q*w^(2^kappa)*delta^q + c^q*w^(2^kappa)*delta^q)*x*y^(2^kappa)+(w^(2^kappa +1)*delta +a*w^(2^kappa + q)*delta + b*w^(2^kappa*q+1)*delta +c*w^(2^kappa*q+q)*delta + w^(2^kappa*q+q)*delta^q + a^q*w^(2^kappa*q+1)*delta^q + b^q*w^(2^kappa+q)*delta^q + c^q *w^(2^kappa+1)*delta^q)*y^(2^kappa +1)"

def returnKim(sol):
	expression = "(1+a+b+c)*x^(2^kappa +1) + (w+a*w^q +b*w +c*w^q)*x^(2^kappa)*y + (w^(2^kappa) + a*w^(2^kappa) + b*w^(2^kappa*q) + c*w^(2^kappa *q))*x*y^(2^kappa) + (w^(2^kappa + 1) +a*w^(2^kappa + q) +b*w^(2^kappa*q +1) +c*w^(2^kappa*q+q))*y^(2^kappa + 1)"
	def Kim((x,y)):
		res = sage_eval(expression, locals={'a':sol[2],'b':sol[3], 'c':sol[4], 'x':x, 'y':y, 'w':w, 'q':q, 'kappa':kappa})
		return res
	return Kim
			
	
def returnFunction(sol, pos):
	'''0 -> gamma, 1 -> beta'''
	def func((x,y)):
		res = sage_eval(genericComponent, locals={'delta':sol[pos],'a':sol[2],'b':sol[3], 'c':sol[4], 'x':x, 'y':y, 'w':w, 'q':q, 'kappa':kappa})
		#print (res in L)
		return res
	return func

cands = []

print("Find subfield linearly independent scalars gamma, beta:")
for beta in tqdm(k):
	for gamma in k:
		for scalar in L:
			if beta == scalar*gamma:
				break
		else:
			cands.append( (gamma,beta) )
print("we have %d such pairs."%len(cands))

def checkSolution(sol):
	kim = returnKim(sol)
	F = returnFunction(sol,0)
	G = returnFunction(sol,1)
	for vec in LtimesL:
		res = K2LtL[ kim(vec) ]
		if (F(vec),G(vec)) != res:
			tqdm.write("Trace function components do not match")
			return False
	else:
		return True

def checkPair(sols, pair):
	cont= True
	gamma = pair[0]
	beta = pair[1]
	for a in tqdm(K):
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
										tqdm.write("First success, found a,b,c...")
										sol = (gamma,beta,a,b,c)
										tqdm.write( str(sol))
										if checkSolution(sol):
											tqdm.write("yeah we got it!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
											sols.append( sol )
											tqdm.write(len(sols))
											cont = False
											break
										
i = 0			
sols = []
while len(sols)<1:
	checkPair(sols,cands[i])
	i+=1

#sol=sols[0]

#sol = (w^2, w, 0, w, w^9 + w^6 + w^4 + w^3 + w^2 + w)



# F_1, F_2 affine mappings as in paper
def returnComp(func,pos):
	def F((x,y)):
		#print ("Handling %s with return comp: %s."%(str((x,y)),(x,y)[pos])) 
		#print func((x,y))
		#print ("\n")
		return ((x,y)[pos], func((x,y)))
	return F


		
#print("output of kim function:")
#for vec in LtimesL:
	#res = kim(vec)
	#print K2LtL[res]
	#print (F(vec), G(vec))
	#print("\n")
	
#F1 = returnComp(F,0)
#F2 = returnComp(G,1)

#F1Table = dict()
#F2Table = dict()


#for vec in tqdm(LtimesL):
	#res = F1(vec)
	#if res not in F1Table.values():F1Table[vec] = res		
#for vec in tqdm(LtimesL):
	#res = F2(vec)
	#if res not in F2Table.values():F2Table[vec] = res
		
#print len(F1Table)
#print len(F2Table)

#F1TableInv = dict()
#for k,v in F1Table.iteritems():
	#F1TableInv[v] = k

#apnMatrix = matrix(K,2,0)
#g = dict()
#for vec in LtimesL:
	##print("%s maps to:"%str(vec))
	#intermediate = F1TableInv[vec]
	#res = F2Table[intermediate]
	
	## construct apn field matrix
	#fieldElem = vec[0] + w*vec[1]
	#resElem = res[0] + w*res[1]
	#M = matrix([[fieldElem],[resElem]])
	#g[fieldElem] = resElem
	#print M
	#apnMatrix = apnMatrix.augment( M )
	##print("%s \n"%str(res))
	
#apnMatrix = apnMatrix[:, 1:]
