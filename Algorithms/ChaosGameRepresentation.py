import math
import numpy as np
import pylab
import sys
from matplotlib import cm
from collections import Counter

coords = {'A': (0, 0), 'C': (0, 1), 'T': (1, 0), 'G': (1, 1)}

def setup(filename):
	with open(filename) as f:
		data = "".join(f.read().split("\n")[1:])
	return data

# returns all kmers and amount of kmers
# deletes kmers with "N"
def count_kmers(sequence, k):
	kmers = []
	toremove = []
	for i in range(len(sequence) - (k - 1)):
		kmers.append(sequence[i: i + k])
	for i in range(len(kmers)):
		if "N" in kmers[i]:
			toremove.append(i)
			print("Found a N uncertainty:", kmers[i])
	for index in reversed(toremove):
		del kmers[index]
	return kmers, len(kmers)

def probability(count, k, length):
	return count / (length - k + 1)

def step(current, aa):
	newX = current[0] + ((1/2) * (coords[aa][0] - current[0]))
	newY = current[1] + ((1/2) * (coords[aa][1] - current[1]))
	return newX, newY

def calcPos(kmer):
	current = (0.5, 0.5)
	for el in kmer:
		current = step(current, el)
	size = math.pow(2, len(kmer))
	current = (current[0] * size, current[1] * size)
	pos = (math.floor(size - current[1]), math.floor(current[0]))
	return pos

def ChaosGameRepresentation(filename, k):
	data = setup(filename)
	length = len(data)
	kmers, size = count_kmers(data, k)
	c = Counter(kmers)
	amounts = c.most_common()
	size = math.floor(math.pow(2, k))
	probabilityMatrix = np.zeros((size, size))
	for am in amounts:
		coords = calcPos(am[0])
		probabilityMatrix[coords] = probability(am[1], k, length)
	pylab.imshow(probabilityMatrix, interpolation='nearest', cmap=cm.gray_r)
	pylab.show()

def main():
	if len(sys.argv) < 3:
		print("Please input your arguments in the following order: \n [script] -f [filename] -k [n kmers]")
		return
	filename = sys.argv[2]
	k = sys.argv[4]
	ChaosGameRepresentation(filename, int(k))

main()
