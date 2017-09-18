from pyeda.inter import *
import pickle

def toCNF(readFile, outFile):
	print("Opening file.")
	with open(readFile, "rb") as f:
		expression1, expressionList = pickle.load(f)
	
	print("Parsing first expression.")
	expression = expr(expression1)
	counter = 1		# for naming auxilliary variables
	length = len(expressionList)
	for e in expressionList:
		print ("Reading expression %d/%d."%(counter, length))
		e = expr(e)
		print("Tseitining.")
		e = e.tseitin("z%d"%counter)
		expression = And( expression , e)
		counter += 1
		
	print("Encoding as CNF.")
	lit, nvars, clauses = expression._encode_cnf()
	
	print("Generating file.")
	final = "p cnf %d %d \n"%(nvars, len(clauses))
	
	# Save cnf
	for clause in clauses:
		final += " ".join(str(idx) for idx in clause) + " 0\n"
	
	print("Saving it under %s."%outFile)
	with open(outFile, "w") as f:
		f.write(final)
	return 
 
		
	
