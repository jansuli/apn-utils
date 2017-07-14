from pyeda.inter import *

def toCNF(readFile, outFile):
	print("Opening file.")
	with open(readFile, "r") as f:
		expression = f.read()
	print("Turning it into pyeda expression.")	
	expression = expr(expression)
	
	print("Applying Tseitin transformation.")
	transformed = expression.tseitin()
	
	print("Getting solver ready string.")
	dimacs = expr2dimacscnf(transformed)
	cnf = (dimacs[1]).__str__()
	#print(dimacs[0])
	
	print("Saving under '%s'."%outFile)
	with open(outFile, "w") as f:
		f.write(cnf)
	print("Done.")
	
	#return transformed, expression
	
