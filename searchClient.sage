from multiprocessing.managers import BaseManager
from multiprocessing import Process,current_process
from Queue import Empty, Full
from time import sleep
m=6
K.<a> = GF(2^m, 'a')
V = K.vector_space()
BigV = VectorSpace(GF(2), m*2)
print("The minimal polynomial of 'a' is:")
print(a.minimal_polynomial())
k = K.list()
upper = matrix(k[1:])

def GFtoBinMatrix(M, dim):
  # iterate over the columns
  A = []
  N = M.ncols()
  for i in range(0, N):
    col = M[:, i]
    tmpCol = []
    for elem in col.list():
      coeff = vector(elem).list()
      tmpCol = tmpCol + coeff
    A.append(tmpCol)
  return matrix(GF(2),A).transpose()

upper = GFtoBinMatrix(upper,m)


def make_client_manager(ip, port, authkey):
    class ServerJobManager(BaseManager):
        pass
    ServerJobManager.register('get_Jobs')
    ServerJobManager.register('get_Solutions')

    manager = ServerJobManager(address=(ip,port), authkey = authkey)
    manager.connect()

    print("Client connected to %s:%s."%(ip,port))
    return manager
    
def check_worker(job_queue, solutions):
    """A worker function to be launched in a separate process. Takes jobs from the job_queue -
    each a tuple of a priority number and another tuple consisting of a matrix constructed so far and options for the next column
    to be constructed."""
    while True:
        try:
            job = job_queue.get_nowait()
            jobTuple = job[1]
            check_options(job_queue, jobTuple, solutions)

        except Empty:
            print("No more numbers to crunch...")
            return

def check_options(job_queue ,jobTuple, solutions):
    mat = jobTuple[0]
    options = jobTuple[1]
    print("------------------------------\n%s\n-------------------------------"%mat.str())
    nCols = mat.ncols()
    if nCols < 3:
        combIndices = list(range(0,nCols+1))
        combIndices = [combIndices]
    else:
		combIndices = Combinations(range(nCols+1), 4)
    combIndices = [ind for ind in combIndices if nCols in ind and len(ind) > 1]

    # generate Lookup with apn-requirement
    testMatrices = list()
    for ind in combIndices:
        cols = mat[:,ind[:-1]]
	testMatrices.append(cols)

    lookUp = list()
    upperNew = upper[:,nCols]
    for option in options:
        lowerNew = matrix([option]).transpose()
        colNew = upperNew.stack(lowerNew)
        for testMatrix in testMatrices:
		    tMat = testMatrix.augment(colNew)
		    r = int(tMat.rank())
		    if r != 4 and nCols > 3:
			    break
		    elif r != tMat.ncols():
			    #print("\nWith \n%s, broke at \n%s, cause rank is %d vs %d"%(mat.str(),tMat.str(), r, nCols+1))
			    break
	else:
	    optionsCopy = list(options)
	    optionsCopy.remove(option)
	    lookUp.append( (colNew, optionsCopy) )

    if nCols < 2^m-2:
        for col, opt in lookUp:
            Mat = mat.augment(col)
            #print("Whilst having %d solutions worker %d appends\n%s\n"%(solutions.qsize(), current_process().pid, Mat.str()))
            try:
			    job_queue.put_nowait( (2^m-2-nCols, (Mat, opt)) )
            except Full:
				print("something went wrong")
    else:
        for col, opt in lookUp:
            Mat = mat.augment(col)
            solutions.put(Mat)
        return "Done"

def initial_jobs(job_queue):
    unused = V[1:]
    for elem in unused:
        lowerCol = matrix([elem]).transpose()
        mat = upper[:,0].stack(lowerCol)
        unusedCopy = list(unused)
        unusedCopy.remove(elem)
        job_queue.put( (2^m-2 , (mat, unusedCopy)) )
        print("putting\n%s\n in queue"%mat)
    sleep(2)
    return

def start_workers(job_queue, solutions, nWorkers = 2):
    workers = []
    for i in range(nWorkers):
        p = Process(target=check_worker, args=(job_queue, solutions))
        workers.append(p)
        p.start()
    for p in workers:
        p.join()

if __name__ == '__main__':
    manager = make_client_manager('localhost', 3333, 'pwd')
    JOBS = manager.get_Jobs()
    print(repr(JOBS))
    SOLS = manager.get_Solutions()
    initial_jobs(JOBS)
    start_workers(JOBS,SOLS, nWorkers = 8)
    
