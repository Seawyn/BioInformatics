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

def setupMatrix(seq1, seq2, gapScore):
	optAlignment = []
	col = [(0, '-')]
	for e in seq1:
		index = len(col)
		col.append(((col[(index - 1)][0] - gapScore), 'l'))
	optAlignment.append(col)
	for e in seq2:
		col = []
		index = len(optAlignment)
		col.append(((optAlignment[(index - 1)][0][0] - gapScore), 't'))
		optAlignment.append(col)
	return optAlignment

def validatePosition(l, t, b, gapScore, ns):
	values = []
	values.append(((l[0] - gapScore), 'l'))
	values.append(((t[0] - gapScore), 't'))
	if (ns[0] > ns[1]):
		values.append(((b[0] + BLOSUM50[indexBLOSUM[ns[0]]][indexBLOSUM[ns[1]]]), 'b'))
	else:
		values.append(((b[0] + BLOSUM50[indexBLOSUM[ns[1]]][indexBLOSUM[ns[0]]]), 'b'))
	maxValue = max(values)
	if (values[2][0] == maxValue[0]):
		return values[2]
	else:
		return maxValue

def traceback(alignmentMatrix):
    trace = []
    i = len(alignmentMatrix) - 1
    j = len(alignmentMatrix[0]) - 1
    cursor = alignmentMatrix[i][j]
    trace = [cursor] + trace
    while (not (i == 0 and j == 0)):
        if (cursor[1] == 'l'):
            j -= 1
        elif (cursor[1] == 't'):
            i -= 1
        else:
            i -= 1
            j -= 1
        cursor = alignmentMatrix[i][j]
        trace = [cursor] + trace
    return trace

def buildAnswer(trace, seq1, seq2):
    seq1 = list(seq1)
    seq2 = list(seq2)
    answer = ["", ""]
    for el in trace[1:]:
        if (el[1] == 'b'):
            answer[0] += seq1.pop(0)
            answer[1] += seq2.pop(0)
        elif (el[1] == 't'):
            answer[0] += '-'
            answer[1] += seq2.pop(0)
        else:
            answer[0] += seq1.pop(0)
            answer[1] += '-'
    print(answer[0])
    print(answer[1])

def Needleman_Wunsch(seq1, seq2, gapScore):
    optAlignment = setupMatrix(seq1, seq2, gapScore)
    for i in range(1, len(seq2) + 1):
        for j in range(1, len(seq1) + 1):
            value = validatePosition(optAlignment[i][j - 1], optAlignment[i - 1][j], 
                    optAlignment[i - 1][j - 1], gapScore, (seq2[i - 1], seq1[j - 1]))
            optAlignment[i].append(value)
    trace = traceback(optAlignment)
    buildAnswer(trace, seq1, seq2)

def main():
	seq1 = "HEAGAWGHEE"
	seq2 = "PAWHEAE"
	gapScore = 8
	populateIndexBLOSUM()
	Needleman_Wunsch(seq1, seq2, gapScore)

main()		
