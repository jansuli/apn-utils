from pyeda.inter import *
import pickle
from tqdm import tqdm

def toCNF(readFile, outFile):
	print("Opening file.")
	with open(readFile, "rb") as f:
		expression1, expressionList = pickle.load(f)
	
	expression2 = expr(expression1)
	counter = 1
	length = len(expressionList)
	for e in expressionList:
		print ("Reading expression %d/%d."%(counter, length))
		e = expr(e)
		print("Tseitining")
		e = e.tseitin("z%d"%counter)
		expression2 = And( expression2 , e)
		counter += 1
		
	print("Encoding")
	lit, nvars, clauses = expression2._encode_cnf()
	
	print("Generating file")
	final = "p cnf %d %d \n"%(nvars, len(clauses))
	for clause in tqdm(clauses):
		final += " ".join(str(idx) for idx in clause) + " 0\n"
	
	with open(outFile, "w") as f:
		f.write(final)
	
	print("Done.")
	return
	#return transformed, expression
	
