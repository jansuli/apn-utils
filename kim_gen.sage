from tqdm import tqdm
from time import sleep
from multiprocessing import Pool, Manager
from numpy import array_split
from functools import partial
import pickle
m = 8
q = 2^(m/2)

K.<w> = GF(2^m, 'w')
L = [K(0)] + [elem for elem in K if elem != 0 and q-1 % elem.multiplicative_order() == 0]	# subfield
k = list(K)
k.remove(K(0))


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
	gamma = pair[0]
	beta = pair[1]
	for a in K:
		#print("checking")
		#print(a)
		for b in K:
			for c in K:
				# Checking x maps to tr^m_k(beta, k(x,y) bijectivity
				# precheck
				if beta + a*beta + b*beta + c*beta in L:
					break
				else:
					left = (beta+a*beta+b*beta+c*beta+beta^q+a^q*beta^q+b^q*beta^q+c^q*beta^q)*(w*beta + a*w^q*beta + b*w*beta +c*w^q*beta + w^q*beta^q+a^q*w^q*beta^q + b^q*w^q*beta^q + c^q*w*beta^q)^(2^kappa)
					right = (beta + a*beta + b*beta + c*beta + beta^q + a^q*beta^q + b^q*beta^q + c^q*beta^q)^(2^kappa) * (w^(2^kappa)*beta +a*w^(2^kappa)*beta+b*w^(2^kappa*q)*beta +c*w^(2^kappa *q)*beta + w^(2^kappa *q)*beta^q + a^q *w^(2^kappa *q)*beta^q + b^q *w^(2^kappa *q)*beta^q + c^q *w^(2^kappa*q)*beta^q)
					
					if left != right:
						#tqdm.write("no success")
						break
						
					else:
						tqdm.write("else clause")
						# check other way round
						# precheck
						if w^(2^kappa + 1)*gamma + a*w^(2^kappa + q)*gamma + b*w^(2^kappa*q+1)*gamma^q + a^q*w^(2^kappa*q+q)*gamma in L:
							break
						else:
							left= (w^(2^kappa +1)*gamma + a*w^(2^kappa + q)*gamma + b*w^(2^kappa*q+1)*gamma + c*w^(2^kappa*q+q)*gamma + w^(2^kappa*q+q)*gamma^q + a^q*w^(2^kappa*q+1)*gamma^q + b^q*w^(2^kappa +q) *gamma^q + c^q * w^(2^kappa+1) *gamma^q)*(w^(2^kappa)*gamma + a*w^(2^kappa)*gamma + b*w^(2^kappa*q)*gamma + c*w^(2^kappa*q)*gamma + w^(2^kappa*q)*gamma^q +a^q *w^(2^kappa*q)*gamma^q + b^q*w^(2^kappa)*gamma^q + c^q*w^(2^kappa)*gamma^q)^(2^kappa)
							right = left
							
							if left != right:
								#tqdm.write("nearly...")
								break
						
			else:
				print("SUCCESS!!!!")
				sols.append((beta,gamma),(a,b,c))
				
if __name__ == "__main__":
	pCount = 8
	chunkCount = len(cands)/(pCount*10)
	chunks = array_split(cands, chunkCount)
	chunks = [arr.tolist() for arr in chunks]
	print chunks[:5]
	m = Manager()
	s = m.list()
	check = partial(checkPair, s)
	p = Pool(processes=pCount)
	for chunk in tqdm(chunks):
		p.map(check,chunk)
		#p.close()
		#p.join()
	p.close()
	
	for sol in s:
		print sol
		
	with open("kim_params.data", "w+") as f:
		pickle.dump(s,f)
				

