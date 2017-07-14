from pyeda.inter import *
import pickle
from tqdm import tqdm

def toCNF(readFile, outFile):
	print("Opening file.")
	newList = []
	with open(readFile, "rb") as f:
		expressionList = pickle.load(f)
	auxCount = 0
	for expression in tqdm(expressionList):
		tqdm.write("\nMaking it pyeda compatible.")
		expression = expr(expression)
		
		if expression.is_cnf():
			tqdm.write("Already cnf. Proceeding.")
			newList+= list(expression.xs)
		else:
			tqdm.write("Applying Tseitin transformation.")
			transformed = expression.tseitin("z%d"%auxCount)
			#print (transformed)
			#print ("\n")
			newList += list(transformed.xs)
			tqdm.write("Done. Proceeding.")
			auxCount += 1
	
	print("Getting solver ready string.")
	whole = expr(And(*newList))
	#wholeExpr = expr(whole)
	dimacs = expr2dimacscnf(whole)
	cnf = (dimacs[1]).__str__()
	#print(dimacs[0])
	
	print("Saving under '%s'."%outFile)
	with open(outFile, "w") as f:
		f.write(cnf)
	print("Done.")
	
	#return transformed, expression
	
