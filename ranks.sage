from tqdm import tqdm
load("classes.sage")

def alternatingBilinearForm(func, vecSpace):
	'''Maps boolean function 'func':F^m -> F (we will later use the component functions) to alternating bilinear form.'''
	
	def Bf(x,y):
		return func(x+y) + func(x) + func(y) + func(vecSpace.zero())
	return Bf
	
def alternatingMatrix(form, basis):
	size = basis[0].parent().dimension()
	mat = matrix(GF(2), size)
	for i in range(1,size):
		for j in range(0,i):
			mat[i,j] = form(basis[i], basis[j])
	mat = mat + mat.transpose()
	return mat
	
test = """x^3,x^9,x^57,g^15*x^48 + g^16*x^33 + g^16*x^18 + x^17 + x^3,x^72 + x^6 + x^3,x^144 + x^6 + x^3,g^21*x^144 + g^183*x^66 + g^245*x^33 + x^3,g^135*x^144 + g^120*x^66 + g^65*x^18 + x^3,g^67*x^192 + g^182*x^132 + g^24*x^6 + x^3,x^160 + x^132 + x^80 + x^68 + x^6 + x^3,x^66 + x^40 + x^18 + x^5 + x^3,x^130 + x^66 + x^40 + x^12 + x^3,g^189*x^192 + g^143*x^144 + g^22*x^132 + g^21*x^129 + g^133*x^96 + g^239*x^72 + g^229*x^66 + g^31*x^48 + g^187*x^36 + g^185*x^33 + g^68*x^24 + g^236*x^18 + g^75*x^12 + g^91*x^9 + g^97*x^6 + g^160*x^3,g^100*x^192 + g^12*x^160 + g^15*x^144 + g^243*x^136 + g^234*x^132 + g^33*x^130 + g^39*x^129 + g^139*x^96 + g^51*x^80 + g^229*x^72 + g^39*x^68 + g^17*x^66 + g^189*x^65 + g^126*x^48 + g^198*x^40 + g^238*x^36 + g^192*x^34 + g^217*x^33 + g^122*x^24 + g^144*x^20 + g^169*x^18 + g^141*x^17 + g^236*x^12 + g^117*x^10 + g^183*x^9 + g^184*x^6 + g^231*x^5 + g^228*x^3,g^155*x^192 + g^96*x^144 + g^223*x^132 + g^77*x^129 + g^88*x^96 + g^232*x^72 + g^69*x^66 + g^142*x^48 + g^168*x^36 + x^33 + g^145*x^24 + g^234*x^18 + g^202*x^12 + g^94*x^9 + g^189*x^6 + g^241*x^3,g^126*x^192 + g^119*x^144 + g^221*x^132 + g^222*x^129 + g^79*x^96 + g^221*x^72 + g^187*x^66 + g^148*x^48 + g^187*x^36 + g^237*x^24 + g^231*x^12 + g^119*x^9 + g^244*x^6 + g^236*x^3,g^151*x^192 + g^13*x^144 + g^58*x^132 + g^143*x^129 + g^110*x^96 + g*x^72 + g^244*x^66 + g^26*x^48 + g^180*x^36 + g^8*x^33 + g^69*x^24 + g^76*x^18 + g^201*x^12 + g^201*x^9 + g^19*x^6 + g^107*x^3,g^86*x^192 + g^224*x^129 + g^163*x^96 + g^102*x^66 + g^129*x^48 + g^102*x^36 + g^170*x^33 + g^14*x^24 + g^170*x^18 + g^101*x^12 + g^58*x^6 + g^254*x^3,g^95*x^192 + g^242*x^144 + g^195*x^132 + g^98*x^129 + g^84*x^96 + g^45*x^72 + g^234*x^66 + g^202*x^48 + g^159*x^36 + g^58*x^33 + g^23*x^24 + g^148*x^18 + g^230*x^12 + g^32*x^9 + g^54*x^6 + g^41*x^3,g^132*x^192 + g^37*x^144 + g^91*x^132 + g^188*x^129 + g^76*x^96 + g^162*x^72 + g^46*x^66 + g^252*x^48 + g^42*x^36 + g^81*x^33 + g^83*x^24 + g^13*x^18 + g^185*x^12 + g^163*x^9 + g^216*x^6 + g^181*x^3,g^91*x^192 + g^124*x^144 + g^214*x^132 + g^106*x^129 + g^59*x^96 + g^172*x^72 + g^138*x^66 + g^163*x^48 + g^58*x^36 + g^100*x^33 + g^32*x^24 + g^250*x^18 + g^45*x^12 + g^241*x^6 + g^157*x^3g^25*x^192 + g^140*x^144 + g^59*x^132 + g^129*x^129 + g^42*x^96 + g^164*x^72 + g^149*x^66 + g^119*x^48 + g^74*x^36 + g^211*x^33 + g^9*x^24 + g^46*x^18 + g^130*x^12 + g^185*x^9 + g^147*x^6 + g^27*x^3,g^113*x^192 + g^56*x^144 + g^68*x^132 + g^155*x^129 + g^91*x^96 + g^78*x^72 + g^159*x^66 + g^30*x^48 + g^194*x^36 + g^14*x^33 + g^238*x^24 + g^91*x^18 + g^100*x^12 + g^96*x^9 + g^222*x^6 + g^178*x^3"""
functions = test.split(',')
dim = 8
K.<a> = GF(2^dim, 'a')
rankDists = dict()
for apn in tqdm(functions):
	f = APN(K,K, apn, 'g')
	f.componentFunctions()
	ranks = dict()
	for i in range(0, dim+1):
		ranks[i] = 0
	tqdm.write("Iterating over component functions...")
	for b,Fb in tqdm(f.components.iteritems()):
		Bf = alternatingBilinearForm(Fb, f.domainVecSpace)
		B = alternatingMatrix(Bf, f.domainVecSpace.basis())
		r = B.rank()
		ranks[r] +=1
		#tqdm.write("%s\n has rank %d."%(B.str(), r))
	tqdm.write("For apn\n%s\n the distrubution is:"%apn)
	tqdm.write(str(ranks))
	rankDists[apn] = ranks
