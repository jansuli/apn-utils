n=6
K.<w> = GF(2^n, 'w', repr="log")
R.<z> = PolynomialRing(K, 'z')
KL = K.list()

#def f(x):
    # Function with non-classical rank distribution on n=8 
    #return w^154*x^192 + w^36*x^160 + w^83*x^144 + w^160*x^136 + w^253*x^132 + w^215*x^130 + w^221*x^129 + w^76*x^96 + w^137*x^80 + w^206*x^72 + w^185*x^68 + w^165*x^66 + w^201*x^65 + w^226*x^48 + w^25*x^40 + w^65*x^36 + w^11*x^33 + w^170*x^24 + w^247*x^20 + w^155*x^18 + w*x^17 + w^146*x^12 + w^204*x^10 + w^121*x^9 + w^202*x^6 + w^246*x^5 + w^170*x^3

def f(x):
    return x^3

def switchingAPN(vec, elem):
    '''Takes a binary vector and field element to return a new APN function according to the Switching Method.'''
    
    def phi(x):
        return f(x) + elem*vec[KL.index(x)]
    return phi

def lagrangeInterpolation(func):
    '''Returns the polynomial corresponding to input "func" via Lagrange-Interpolation.'''
    
    poly = R(0)
    for i in range(2^n):
        ai = KL[i]
        bi = func(ai)
        prod = R(1)
        for k in range(2^n):
            if k != i:
                ak = KL[k]
                prod = prod* ((z-ak)/(ai-ak))
        poly = poly + bi*prod
    return poly
        
switchingAPNs = []
for g in KL:
    if g!= 0:
        M = matrix(GF(2), 0, 2^n)
        rank = 0
        print("g=%s. Searching constraints."%str(g))
        for x in K:
            for y in K:
                for a in K:
                    if f(x) + f(x+a) + f(y) + f(y+a) == g:
                        constraintTuple = (x, x+a, y, y+a)
                        newRow = matrix(GF(2), 1, 2^n)

                        for elem in constraintTuple:
                            newRow[0,KL.index(elem)] = 1
                        
                        if not vector(newRow) in M.rows():
                            M = M.stack(newRow)
                            rank = M.rank()
                            print("Added row to M which now has rank %d."%rank)
                            
                            if rank >= 2^n:
                                print("Rank of M is too big, only solution is zero.")
                                break
                else:
                    continue
                break
            else:
                continue
            break
        else:
        # Successfully built constraint matrix M for given g.
            MKernel = M.right_kernel()
            print("We have %d Boolean functions."%(2^MKernel.dimension()-1))
            
            #print("Now interpolating them.")
            #newAPNs = [lagrangeInterpolation(switchingAPN(vec,g)) for vec in M.right_kernel()]
            #switchingAPNs += newAPNs                              
        break
