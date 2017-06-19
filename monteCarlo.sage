from time import time
import datetime
from functools import partial
import multiprocessing as mp

def checkMove(mat, nCols,V, w, inds,move):
	newColUp = matrix(GF(2),vector(w^nCols)).transpose()
	newColLo = matrix(GF(2),vector(move)).transpose()
	newCol = newColUp.stack(newColLo)
	matNew = mat.augment(newCol)
		
	N = nCols + 1
	if N <= 3:
		cols = matNew.columns()
		if not V.are_linearly_dependent(cols):
			return move
	else:
		#print matNew.str()
		for ind in inds:
			cols = matNew[:, ind].columns()
			if V.are_linearly_dependent(cols):	# faster than rank checking it seems
				return
		else:
			return move


class Board():
	def __init__(self, m = 3, nWorkers = 2):
		self.state = ()
		K.<w> = GF(2^m, 'w')
		V = VectorSpace(GF(2), 2*m)
		self.nWorkers = nWorkers
		
		self.field = K
		self.vecSpace = V
		
		self.legalDict = dict()
		self.combIndices = dict()
		for i in range(4,2^m):
			self.combIndices[i] = [ind for ind in Combinations(i, 4).list() if i-1 in ind]
		
	def start(self):
		# Return starting state, i.e. empty tuple
		return ()
		
	def make_move(self, move):
		state = list(self.state)
		state.append(move)
		self.state = tuple(state)
	
	def next_state(self, state, move):
		# Takes game state and move to be applied
		current = list(state)
		current.append(move)
		return tuple(current)
	
	def legal_plays(self, state_history):
		current = state_history[-1]
		
		K = self.field
		if current in self.legalDict:
			print("We know that state...")
			legal = self.legalDict[current]
			return legal
			
		if len(current) == 0:
			legal = K.list()[1:]
			self.legalDict[current] = legal
			return legal
						
		m = K.degree()
		w = K.primitive_element()
		V = self.vecSpace
			
		upper = []
		i = 0
		for elem in current:
			upper.append(w^i)
			i += 1
		
		currentList = list(current)
		upperMatGF = matrix(K, upper)
		lowerMatGF = matrix(K,currentList)
		matGF = upperMatGF.stack(lowerMatGF)
		mat = matrix(GF(2), 2*m, 0)
		
		for col in matGF.columns():
			colUp = matrix(vector(col[0])).transpose()
			colLo = matrix(vector(col[1])).transpose()
			newCol = colUp.stack(colLo)
			mat = mat.augment(newCol)
		nCols = matGF.ncols()
				
		options = [elem for elem in K if elem not in currentList and elem != K(0)]
		legal = []
		
		print("Checking to %d options."%len(options))# \n%s\n"%(len(options),mat.str()))
		
		#for move in options:
			#newColUp = matrix(GF(2),vector(w^nCols)).transpose()
			#newColLo = matrix(GF(2),vector(move)).transpose()
			#newCol = newColUp.stack(newColLo)
			#matNew = mat.augment(newCol)
			
			#N = nCols + 1
			#if N <= 3:
				#cols = matNew.columns()
				#if not V.are_linearly_dependent(cols):
					#legal.append(move)
			#else:
				#inds = self.combIndices[N]
				##print matNew.str()
				#for ind in inds:
					#cols = matNew[:, ind].columns()
					#if V.are_linearly_dependent(cols):	# faster than rank checking it seems
						#break
				#else:
					#legal.append(move)
		if __name__ == "__main__": 
			inds = None if nCols < 3 else self.combIndices[nCols + 1]
			moveChecker = partial(checkMove, mat, nCols, V, w, inds)
			p = mp.Pool(processes = self.nWorkers)
			res = p.map(moveChecker, options)
			legal = [opt for opt in res if opt != None]
			p.close()
			p.join()
			
		self.legalDict[current] = legal
		
		if len(legal) > 0:
			return legal
		else:
			return False		
		
	def win(self, state_history):
		# Takes entire game history and returns 1 if game is won or 0 if game is ongoing.
		m = self.field.degree()
		maxCols = 2^m - 1
		current = state_history[-1]
		if len(current) == maxCols:
			print ("Won!!!!!!!!!!!!!%d"%(2^m -1))
			print current
			with open("monte_perm%.1f.tuple"%time(), "w") as f:
				f.write(str(current))
			return 1
		else:
			return len(current)/maxCols
		
class MonteCarlo(object):
	
	def __init__(self, board, **kwargs):
		self.board = board
		self.states = []
		
		secs = int(kwargs.get("time", 60))
		self.duration = datetime.timedelta(seconds = secs)
		self.maxCols = kwargs.get("maxCols")
		self.C = kwargs.get("C", 1.4)
		
		self.wins = {}
		self.plays = {}
				
	def update(self, state):
		self.states.append(state)
		
	def get_play(self):
		self.max_depth = 0
		state = self.states[-1]
		legal = self.board.legal_plays(self.states[:])
		
		if not legal:
			return False
		if len(legal) == 1:
			return legal[0]
			
		games = 0
		begin = datetime.datetime.utcnow()
		while datetime.datetime.utcnow()-begin < self.duration:
			self.run_simulation()
			games += 1
			
		moves_states = [ (p, self.board.next_state(state, p)) for p in legal ]
		
		print ("Simulated %d games."%games)
		
		percent_wins, move = max(
			(self.wins.get(S, 0) /
			self.plays.get(S,1), p)
			for p,S in moves_states
		)
		
		print "Maximum depth searched: %d."%self.max_depth
		#print self.plays, self.wins
		
		return move
		
	def run_simulation(self):
	
		plays, wins = self.plays, self.wins
		
		visitedStates = set()
		 
		statesCopy = list(self.states)
		state = statesCopy[-1]
		maxMoves = self.maxCols - len(state)
		
		expand = True
		for t in range(1,maxMoves+1):	
			legal = self.board.legal_plays(statesCopy)
			if legal and len(legal)>0:
				moves_states = [ (p, self.board.next_state(state, p)) for p in legal ]
				
				if all(plays.get(S) for p,S in moves_states):
					# if we have stats on all legal moves, use them
					log_total = numerical_approx(log(sum(plays[S] for p,S in moves_states)))
					value, move, state = max(
						(numerical_approx((wins[S] / plays[S]) + self.C * sqrt(log_total / plays[S])), p,S)
						for p,S in moves_states
					)
				else:
					move, state = choice(moves_states)
				#print move, state
					
				statesCopy.append(state)
				
				if expand and state not in self.plays:
					expand = False
					self.plays[state] = 0
					self.wins[state] = 0
					if t > self.max_depth:
						self.max_depth = t
				
				visitedStates.add(state)
				
				win = self.board.win(statesCopy)
				
				if win:
					break
			else:
				break
				
		for state in visitedStates:
			if state not in self.plays:
				continue
			self.plays[state] += 1
			if win:
				self.wins[state] += 1
m = 5
game = Board(m, nWorkers = mp.cpu_count())
monte = MonteCarlo(game, maxCols = 2^m -1, time = 120)
monte.update(game.start())

won = False
while not won:
	move = monte.get_play()
	print("Appending:")
	print(move)
	if move and len(game.state) < 2^m-2:
		game.make_move(move)
		monte.update(game.state)
	elif move and len(game.state) == 2^m - 2:
		game.make_move(move)
		monte.update(game.state)
		print("Won the game :) ")
		won = True
	else:
		# Lost, reset
		print("We lost :(")
		monte.update(game.start())
		game.state = ()
print(game.state)
		
		
