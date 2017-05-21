from multiprocessing.managers import SyncManager
from thread import start_new_thread
from Queue import PriorityQueue,Queue
from time import sleep,time
import pickle
from sage.all import *

def start_Server_JobManager(port):
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
            print(sols)
            with open("apn_perms_%.2f.data"%time(),"w+") as f:
                pickle.dump(sols,f)
        sleep(10)
print_Server_Stats(3333)
