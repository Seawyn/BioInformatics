import math
import numpy as np

trsMatrix = [[0.6, 0.4, 0], [0.25, 0.5, 0.25], [0.25, 0.25, 0.5]]
emProb = [[0.4, 0, 0.3, 0.3], [0.1, 0.4, 0.4, 0.1], [0.4, 0.3, 0, 0.3]]
correspondence = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

S = "CATGCGGGTTATAAC"

path = np.zeros((len(trsMatrix), len(S))).tolist()
forwardProbs = np.zeros((len(trsMatrix), len(S))).tolist()

#returns Es(Emission)
def fetchEmission(state, emission):
	return emProb[(state - 1)][correspondence[emission]]

#for debugging purposes
def printArray(array):
	for line in array:
		print(line)

#log function
def logValues(n):
	if n == 0:
		return None
	else:
		return math.log(n)

#probability of a function
def updateForward(coords):
	prob = fetchEmission(coords[0] + 1, S[coords[1]])
	if prob == 0:
		forwardProbs[coords[0]][coords[1]] = None
	else:
		prob = logValues(prob)
		summation = 0
		for i in range(len(emProb)):
			fp = forwardProbs[i][coords[1] - 1]
			ap = trsMatrix[i][coords[0]]
			if (fp != None):
				#undo ln for sum
				summation += ap * math.pow(math.e, fp)
		if (summation == 0):
			forwardProbs[coords[0]][coords[1]] = None
		else:
			total = prob + math.log(summation)
			forwardProbs[coords[0]][coords[1]] = total

#sum of last calculations 
def Forward():
	prob = 0
	for el in range(len(emProb)):
		if (forwardProbs[el][len(S) - 1] != None):
			prob += math.pow(math.e, forwardProbs[el][len(S) - 1])
	print("Probability of Sequence:", prob)

def Viterbi(S, trsMatrix, emProb):

	#for each column
	for i in range(len(S)):
		if i == 0:

			#for each line in column
			for j in range(len(emProb)):
				path[j][i] = (logValues(fetchEmission((j + 1), S[i]) * 1 / len(emProb)), -1)
				forwardProbs[j][i] = logValues(fetchEmission((j + 1), S[i]) * 1 / len(emProb))
		else:

			#for each line in column
			for j in range(len(emProb)):
				m = (None, 0)
				updateForward((j, i))

				#for each value in max
				for n in range(len(emProb)):
					if path[n][i - 1][0] == None:
						continue
					update = fetchEmission((j + 1), S[i]) * trsMatrix[n][j]
					if update == 0:
						continue
					current = path[n][i - 1][0] + logValues(update)
					if m[0] == None:
						m = (current, n)
					elif m[0] != None and m[0] < current:
						m = (current, n)
				path[j][i] = m

	#calculate path
	states = []
	indexes = []
	for c in range(len(trsMatrix)):
		indexes.append(path[c][len(S) - 1])
	maxInd = indexes[0]
	for ind in indexes:
		if maxInd[0] == None:
			maxInd = ind
		elif ind[0] == None:
			continue
		elif ind[0] > maxInd[0]:
			maxInd = ind
	steps = [indexes.index(maxInd)]
	current = maxInd
	for i in range(len(S) - 2, -1, -1):
		steps = [maxInd[1]] + steps
		if maxInd[1] != -1:
			maxInd = path[maxInd[1]][i]

	#correct offset
	for i in range(len(steps)):
		steps[i] += 1
	print("Most likely sequence of states:", steps)

def main():
	Viterbi(S, trsMatrix, emProb)
	Forward()

main()
