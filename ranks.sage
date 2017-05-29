from multiprocessing import Process, Manager
from Queue import Empty
from time import sleep, time, strftime
import pickle
from jinja2 import Environment, FileSystemLoader
import paramiko
from scp import SCPClient
from funcs import functions
from imtCredentials import *

env = Environment(loader=FileSystemLoader(r"./"))

template = env.get_template('template.html')

load("classes.sage")
load("functions_.sage")

nWorkers = 6
	
dim = 8
K.<a> = GF(2^dim, 'a')

def calcRankDist(jobs, rankQueue,m):
	while True:
		try:
			apn = jobs.get_nowait()
			f = APN(K,K, apn, 'g')
			f.componentFunctions()
			ranks = dict()
			for i in range(0, dim+1):
				ranks[i] = 0
			print("Iterating over component functions...")
			for b,Fb in f.components.iteritems():
				Bf = alternatingBilinearForm(Fb, f.domainVecSpace)
				B = alternatingMatrix(Bf, f.domainVecSpace.basis())
				r = B.rank()
				ranks[r] +=1
				#print("%s\n has rank %d."%(B.str(), r))
			ranks = rankDistCleanUp(ranks)
			print("For apn\n%s\nthe distrubution is:"%apn)
			print(str(ranks))
			
			codeDist = str(codeDistanceFromRanks(ranks, m))
			Walsh = str(WalshFromRanks(ranks,m))
			rankQueue.put( (apn, ranks, codeDist) )
		except SyntaxError:
			rankQueue.put( (apn,"Threw an error...", "") )
		except Empty:
			print("nothing to do")
			break

def updateHomepage(distributions):
	while True:
		results= []
		dists = {}
		for i in range(distributions.qsize()):
			result = distributions.get()
			distribution = str(result[1])
			res = {'apn':result[0], 'distribution':distribution, 'codeDist' : result[2]}#, 'Walsh' : result[3]}
			results.append(res)
			if distribution in dists.keys():
				dists[distribution] += 1
			else:
				dists[distribution] = 1
			
			distributions.put(result)
		t = time()
		upTime = strftime("%d.%m.%Y, %H:%M:%S")
		with open("index.html", "w+") as f:
			f.write(template.render(updateTime =upTime, results = results, distN = dists, cols = len(results)))
		
		ssh = paramiko.SSHClient()
		ssh.load_system_host_keys()
		ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
		ssh.connect('sshgate.uni-paderborn.de', username=imtUser, password=imtPwd)

		with SCPClient(ssh.get_transport()) as scp:
			scp.put('index.html', '~/public/http/index.html')
		sleep(60)

if __name__ == "__main__":
	manager = Manager()
	rankDists = manager.Queue()
	jobs = manager.Queue()
	functions = functions.replace("\n","").split(',')
	for fun in functions:
		jobs.put(fun)
	workers = []
	print(jobs.qsize())
	
	service = Process(target=updateHomepage, args=(rankDists,))
	service.start()
	
	for i in range(nWorkers):
		p = Process(target=calcRankDist, args=(jobs, rankDists, dim))
		p.start()
		workers.append(p)
	for p in workers:
		p.join()
		
	print("Done with work")
	sleep(60)
	service.terminate()
	
	with open("rankQueue%.2f.data"%time(),"w+") as f:
		pickle.dump(rankDists, f)
