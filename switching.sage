n=8
K.<w> = GF(2^n, 'w', repr="log")
R.<z> = PolynomialRing(K, 'z')
KL = K.list()

def f(x):
    return w^154*x^192 + w^36*x^160 + w^83*x^144 + w^160*x^136 + w^253*x^132 + w^215*x^130 + w^221*x^129 + w^76*x^96 + w^137*x^80 + w^206*x^72 + w^185*x^68 + w^165*x^66 + w^201*x^65 + w^226*x^48 + w^25*x^40 + w^65*x^36 + w^11*x^33 + w^170*x^24 + w^247*x^20 + w^155*x^18 + w*x^17 + w^146*x^12 + w^204*x^10 + w^121*x^9 + w^202*x^6 + w^246*x^5 + w^170*x^3

def switchingAPN(vec, elem):
    def phi(x):
        return f(x) + elem*vec[KL.index(x)]
    return phi

def lagrangeInterpolation(func):
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

def algebraicDegree(pol):
    terms = str(pol).replace(" ", "").split("+")
    degrees = []
    for term in terms:
        subTerms = term.split("*")
        if len(subTerms) > 1:
            monom = subTerms[1]
        else:
            if 'z' in term:
                monom = term
            else:
                monom = "z^0"
        exp = 1 if "^" not in monom else int(monom.split("^")[1])
        deg = bin(exp).count('1')

        degrees.append(deg)
        
    return max(degrees)
                              
def buildConstraintMatrix(constraintIndexSet):
    M = matrix(GF(2), len(constraintIndexSet), 2^n)
    row = 0
    for constraintIndexTuple in constraintIndexSet:
        positions = []
        for elem in constraintIndexTuple:
            positions.append(KL.index(elem))
        M[row, positions] = 1
        row += 1
    return M
        
switchingAPNs = []
for g in KL:
    if g!= 0:
        constraintTuples = set()
        print("g=%s. Searching constraints."%str(g))
        for x in K:
            for y in K:
                for a in K:
                    if f(x) + f(x+a) + f(y) + f(y+a) == g:
                        constraintTuples.add( Set([x,x+a, y, y+a]))
        M = buildConstraintMatrix(constraintTuples)
        print("Looking for switching APNs with g=%s."%str(g))
        newAPNs = [lagrangeInterpolation(switchingAPN(vec,g)) for vec in M.right_kernel()]
        print("Found %d new APNs.")
        switchingAPNs += newAPNs                              
