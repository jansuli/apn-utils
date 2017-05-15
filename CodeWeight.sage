from sage.coding.binary_code import weight_dist

def gold(fieldElem):
  return fieldElem**(9)

def coeffList(dim, fieldElem):
  p = fieldElem.polynomial().list()
  padding = [0 for j in range(0, dim - len(p)) ]

  return p+padding

def buildGenMatrix(K, dim):
  M = []
  for e in K:
    res = gold(e)

    coeff1 = coeffList(dim, e)
    coeff2 = coeffList(dim, res)

    M.append ( [1] + coeff1 + coeff2 )
  return matrix(GF(2), M).transpose()

