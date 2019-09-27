import copy

BLOSUM50 = [
[5],
[-1, 13],
[-2, -4, 8],
[-1, -3, 2, 6],
[-3, -2, -5, -3, 8],
[0, -3, -1, -3, -4, 8],
[-2, -3, -1, 0, -1, -2, 10],
[-1, -2, -4, -4, 0, -4, -4, 5],
[-1, -3, -1, 1, -4, -2, 0, -3, 6],
[-2, -2, -4, -3, 1, -4, -3, 2, -3, 5],
[-1, -2, -4, -2, 0, -3, -1, 2, -2, 3, 7],
[-1, -2, 2, 0, -4, 0, 1, -3, 0, -4, -2, 7],
[-1, -4, -1, -1, -4, -2, -2, -3, -1, -4, -3, -2, 10],
[-1, -3, 0, 2, -4, -2, 1, -3, 2, -2, 0, 0, -1, 7],
[-2, -4, -2, 0, -3, -3, 0, -4, 3, -3, -2, -1, -3, 1, 7],
[1, -1, 0, -1, -1, 0, -1, -3, 0, -3, -2, 1, -1, 0, -1, 5],
[0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 2, 5],
[0, -1, -4, -3, -1, -4, -4, 4, -3, 1, 1, -3, -3, -3, -3, -2, 0, 5],
[-3, -5, -5, -3, 1, -3, -3, -3, -3, -2, -1, -4, -4, -1, -3, -4, -3, -3, 15],
[-2, -3, -3, -2, 4, -3, 2, -1, -2, -1, 0, -2, -3, -1, -1, -2, -2, -1, 2, 8]]

indexBLOSUM = {}

def populateIndexBLOSUM():
	population = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R',
					'S', 'T', 'V', 'W', 'Y']
	for p in population:
		indexBLOSUM[p] = len(indexBLOSUM)

def setupMatrix(seq1, seq2):
    optAlignment = []
    col = [(0, '-')]
    for e in seq1:
        col.append((0, 'l'))
    optAlignment.append(col)
    for e in seq2:
        col = [(0, 't')]
        optAlignment.append(col)
    return optAlignment

def validatePosition(pvScores, gapScore, ns):
	values = []
	order = ''
	value = (0, '-')
	values.append((pvScores[0] - gapScore, 't'))
	values.append((pvScores[1] - gapScore, 'l'))
	if (ns[0] > ns[1]):
		valueBLOSUM = BLOSUM50[indexBLOSUM[ns[0]]][indexBLOSUM[ns[1]]]
	else:
		valueBLOSUM = BLOSUM50[indexBLOSUM[ns[1]]][indexBLOSUM[ns[0]]]
	values.append((pvScores[2] + valueBLOSUM, 'b'))
	maxValue = max(values)
	if (maxValue[0] < 0):
		return value
	else:
		for el in values:
			if (el[0] == maxValue[0]):
				order += el[1]
		return (maxValue[0], order)

def traceback(matrix, ii, ij):
	i = ii
	j = ij
	optAlignment = copy.deepcopy(matrix)
	options = []
	trace = []
	cursor = optAlignment[i][j]
	while (cursor[0] != 0):
		trace = [cursor] + trace
		if (cursor[1] == 'l'):
			j -= 1
		elif (cursor[1] == 't'):
			i -= 1
		elif (cursor[1] == 'b'):
			i -= 1
			j -= 1
		elif (len(cursor[1]) > 1):
			for e in cursor[1]:
				optAlignment[i][j] = (optAlignment[i][j][0], e)
				options.append(traceback(optAlignment, ii, ij))
			break
		cursor = optAlignment[i][j]
	if (options == []):
		return trace
	else:
		return options

def buildAnswer(trace, i, j, seq1, seq2):
	algn1 = ""
	algn2 = ""
	for e in range(len(trace)):
		cursor = trace.pop()
		if (cursor[1] == 'l'):
			algn1 = seq1[j - 1] + algn1
			algn2 = '-' + algn2
			j -= 1
		elif (cursor[1] == 't'):
			algn1 = '-' + algn1
			algn2 = seq2[i - 1] + algn2
			i -= 1
		elif (cursor[1] == 'b'):
			algn1 = seq1[j - 1] + algn1
			algn2 = seq2[i - 1] + algn2
			i -= 1
			j -= 1
	print(algn1)
	print(algn2)

def findElement(matrix, el):
	for i in range(len(matrix)):
		for j in range(len(matrix[0])):
			if (matrix[i][j] == el):
				return (i, j)
	return -1

def formatTraces(options):
	traces = []
	for el in options:
		if (isinstance(el[0], list)):
			for j in el:
				traces.append(j)
		else:
			traces.append(el)
	return traces

def Smith_Waterman(seq1, seq2, gapScore):
    optAlignment = setupMatrix(seq1, seq2)
    maxValue = 0
    for i in range(1, len(seq2) + 1):
        for j in range(1, len(seq1) + 1):
            pvScores = [optAlignment[i - 1][j][0], optAlignment[i][j - 1][0], optAlignment[i - 1][j - 1][0]]
            score = validatePosition(pvScores, gapScore, (seq1[j - 1], seq2[i - 1]))
            optAlignment[i].append(score)
            if (score[0] > maxValue):
                maxValue = score[0]

    maxCoords = []
    for i in range(len(optAlignment)):
        for j in range(len(optAlignment[0])):
            if (optAlignment[i][j][0] == maxValue):
                maxCoords.append((i, j))
    for coords in maxCoords:
        traces = traceback(optAlignment, coords[0], coords[1])
        traces = formatTraces(traces)
        if (isinstance(traces[0], list)):
            for el in traces:
                subcoords = findElement(optAlignment, el[len(el) - 1])
                buildAnswer(el, subcoords[0], subcoords[1], seq1, seq2)
        else:
            buildAnswer(traces, coords[0], coords[1], seq1, seq2)

def main():
	seq1 = "WPIWPC"
	seq2 = "IIWPI"
	gapScore = 8
	populateIndexBLOSUM()
	Smith_Waterman(seq1, seq2, gapScore)

main()
