from tqdm import tqdm
m = 6
q = 2^(m/2)

load("functions_.sage")

K.<w> = GF(2^m, 'w')

def kim(params):
	(a,b,c) = params
	def Kim(x):
		return x^3 + a*x^(2+q) + b*x^(2*q+1) + c*x^(3*q)
	return Kim
	
verify = kim((K(1), K(1), w))
print checkApnByDefinition(verify, K)
	
for a in tqdm(K):
	for b in tqdm(K):
		for c in tqdm(K):
			Kim = kim( (a,b,c) )
			if checkApnByDefinition(Kim, K):
				tqdm.write("Params that lead to an APN: %s"%(str((a,b,c))))
	