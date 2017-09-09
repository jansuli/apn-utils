dim = 8
k = 4
c = 2^(7) 	# c as used by Edel, adapted for n=8, k=2

def generalizedKrawtchoukPol(j, i):
	'''Evaluates the j-th generalized Krawtchouk-Polynomial in i.''' 
	result = 0
	for l in range(j+1):
		result += (-1)^(j-l) * 4^(binomial(j-l,2)) * bnaryGaussianCoeff(k-l, k-j) * bnaryGaussianCoeff(k-i, l) *c^l
	return result
	
def bnaryGaussianCoeff(x,k, b=4):
	'''Returns the b-nary Gaussion Coefficient for 'x choose k'.'''
	if k != 0:
		result = 1
		for i in range(k):
			result = result * ((b^x-b^i)/(b^k-b^i))
		return result
	else:
		return 1

# Find all possible rank distributions.
distributions = []
a_0 = 1
for a_1 in range(6):
    for a_2 in range(18):
        for a_3 in range(86):
            if 21*a_1 + 5*a_2 + a_3 == 85:
                a_4 = 255-(a_1+a_2+a_3)
                if 170 <= a_4 <= 240:
                    distributions.append( (a_0, a_1, a_2, a_3, a_4) )
                    
	
print("There are %d possible rank distributions."%len(distributions))
	
# Generate Transformation-Matrix P
P = matrix(ZZ, k+1, k+1)
for i in range(k+1):
	for j in range(k+1):
		P[i,j] = generalizedKrawtchoukPol(j,i)

# Calculate dual distributions.		
dualRanks = dict()
for rankDist in distributions:
	distVector = vector(ZZ, rankDist)
	dualRanks[rankDist] =  2^(-dim) * distVector *P 
	
# Print everything
print("For n=8, the distributions (a_0, …, a_4) and their dual distributions (a_0', …, a_4') are:")
for dist, dual in dualRanks.iteritems():
	print("%20s and its dual distribution %20s."%(str(dist), str(dual)))
