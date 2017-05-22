from multiprocessing.managers import SyncManager
from thread import start_new_thread
from Queue import PriorityQueue,Queue
from time import sleep,time
import pickle
from sage.all import *

def start_Server_JobManager(port, getSaved = False):
    job_queue = PriorityQueue()
    solutions = Queue()
    class JobManager(SyncManager):
        pass
    JobManager.register('get_Jobs', callable=lambda: job_queue)
    JobManager.register('get_Solutions',callable=lambda: solutions)
    manager = JobManager(address=('localhost',port), authkey = "pwd")
    #print("test")
    s = manager.get_server()
    s.serve_forever()
    
manager = start_new_thread(start_Server_JobManager, (3333,))
def print_Server_Stats(port):
    try:
		class JobManager(SyncManager):
			pass

		JobManager.register('get_Jobs')
		JobManager.register('get_Solutions')
		manager = JobManager(address=('localhost',port ), authkey = "pwd")
		manager.connect()
		solutions = manager.get_Solutions()
		job_queue = manager.get_Jobs()
		while True:
			nSols = solutions.qsize()
			nJobs = job_queue.qsize()
			print("At the moment we got %d solutions & %d jobs."%(nSols,nJobs))
			if nSols > 0:
				sols = list()
				for i in range(nSols):
					sol = solutions.get()
					sols.append(sol)
					solutions.put(sol)
				with open("apn_perms_%.2f.data"%time(),"w+") as f:
					pickle.dump(sols,f)
			sleep(2)
    except KeyboardInterrupt:
        print("Saving jobs...")
        jobs = list()
        i=1
        while job_queue.qsize()>0:
			job = job_queue.get_nowait()
			jobs.append(job)
			print("saved %d jobs."%i)
			i+=1
        with open("saved_jobs.tmp","w+") as f:
		    pickle.dump(jobs, f)
        print("saving solutions...")
        sols = list()
        while solutions.qsize()>0:
		    sol = solutions.get_nowait()
		    sols.append(sol)
        with open("saved_sols.tmp","w+") as f:
		    pickle.dump(sols, f)
        return	

print_Server_Stats(3333)
