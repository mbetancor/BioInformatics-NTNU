"""
Documentations:

* Illumina Adapter Sequences: http://dnatech.genomecenter.ucdavis.edu/wp-content/uploads/2013/06/illumina-adapter-sequences_1000000002694-00.pdf
* ABI Ion S5: 

More info (BioInformatics):
https://bioinformatics.amc.nl/wp-content/uploads/gs-sequence-analysis/2018-sequencing-techniques.pdf

"""


import sys, zipfile, gzip, math, random
import matplotlib.pyplot as plt
import numpy as np
from suffixTree import *

def main():
	a = "TGGAATTCTCGGGTGCCAAGGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"
	task = sys.argv[1]
	file = sys.argv[2]
	ans = ""

	if (task=="task1"):
		adapter = a[::-1]
		ans = task1(file, adapter)
	elif (task=="task2"):
		print("hei")
		adapter = a
		ans = task2(file, adapter)
	elif (task=="task22"):
		print("hei")
		adapter = a
		ans = task22(file, adapter)
	elif (task=="task3"):
		print("prueba")
		ans = task3(file)
	elif (task=="task31"):
		ans = task31(file)
	elif(task=="task4"):
		print("prueba4")
		ans=task4(file)
	elif(task=="task44"):
		print("prueba4b")
		ans=task44(file)

	return ans 

def task1(file, adapter):
	
	adapter_tree = SuffixTree()
	adapter_tree.updateTree(adapter)
	hits = 0
	remainingLengths = dict([])
	if (file=="s_3_sequence_1M.txt"):
		lineLen = 50
	elif (file=="s_1-1_1M.txt"):
		lineLen = 76
	elif (file=="seqset3.txt"):
		lineLen = 51

	with open(file, 'r') as f:
		for line in f:
			line = line.rstrip()
			line = line[::-1]
			stuffInTheTree = adapter_tree.isInTree(line) 
			if (stuffInTheTree != []):
				hits += 1
				suffix = stuffInTheTree[-1]
				remainingLength = lineLen - len(suffix)
				if remainingLength not in remainingLengths:
					remainingLengths[remainingLength] = 0
				else:
					remainingLengths[remainingLength] += 1

	del adapter_tree
	print(remainingLengths)
	print(hits) # should output 646364

	
	return hits

def task2(file, adapter):
	
	hits = 0
	counter = 0
	adapter = adapter
	remainingLengths = dict([])
	with open(file, 'r') as f:
		for line in f:

			line = line.rstrip()
			#line = line[-maxLen:]
			r = readSequence(line,adapter, 0.9)
			if (r):
				hits+=1
				remainingLength = 50 - len(r)
				if remainingLength not in remainingLengths:
					remainingLengths[remainingLength] = 0
				else:
					remainingLengths[remainingLength] += 1
			counter=counter+1
			#print(hits)
			print("progress is: " + str(counter/1000000))
			# guardar este valor para luego
	print(hits) # should output 646364
	print(remainingLengths)
	
	return hits

def task22(file, adapter):
	
	# adapter_tree = SuffixTree()
	# adapter_tree.updateTree(adapter)

	hits = 0
	counter = 0
	adapter = adapter[:50]
	remainingLengths = dict([])
	with open(file, 'r') as f:
		for line in f:

			line = line.rstrip()
			#line = line[-maxLen:]
			r = readSequence(line,adapter, 0.75)
			if (r):
				hits+=1
				remainingLength = 50 - len(r)
				if remainingLength not in remainingLengths:
					remainingLengths[remainingLength] = 0
				else:
					remainingLengths[remainingLength] += 1
			counter=counter+1
			#print(hits)
			print("progress is: " + str(counter/1000000))
			# guardar este valor para luego
	print(hits) # should output 646364
	print(remainingLengths)
	
	return hits

def task4(file):
	adapter = "TCGTATGCCGTCTTCTGCTTGAAA"
	barcodes = dict([])
	with open(file,'r') as f:
		for line in f:
			line = line.rstrip()
			#line.replace(adapter,"")
			auxLine = line[-24:]
			line = line[:-24]
			if auxLine==adapter:
				for index in range(4,9):
					if line[-index:] not in barcodes:
						barcodes[line[-index:]]=0
					else:
						barcodes[line[-index:]]+=1
	barcodes = {sequence: occurrences for sequence, occurrences in barcodes.items() if occurrences > 50000}
	print(barcodes)

def task44(file):
	adapter 	= "TCGTATGCCGTCTTCTGCTTGAAA"
	barcodes 	= ["TAGGA", "TAACG", "TATCA", "TCATC", "TCACT", "TCTCC", "TTACA", "TCCGT", "TCCTA", "TCGCG", "TGATG", "TCCAC"]
	sequencesForBarcode = dict([])
	for barcode in barcodes:
		sequencesForBarcode[barcode]=[]
	sequences = dict([])
	with open(file,'r') as f:
		for line in f:
			line = line.rstrip()
			auxLine= line[-24:]
			line = line[:-24]
			barcode = line[-5:]
			sequence = line[:-5]
			if auxLine==adapter and barcode in barcodes:
				if sequence not in sequences:
					sequences[sequence] = 0
				else:
					sequences[sequence] += 1
	sequences = {seq:occ for seq, occ in sequences.items() if occ>100000}
	print(sequences)

def oldtask31(file):
	
	pRandomUnits = 1000000
	indexesOfRandomSamples = random.sample(range(1000000),pRandomUnits)
	minimumOfOcurrences = 10000
	hits = 0
	counter = 0
	sequenceOcurrences = dict([])
	interestingLengths = range(68,77)

	f = open(file)
	lines = f.readlines()
	for index in indexesOfRandomSamples:
		line = lines[index].rstrip()
		for length in interestingLengths:
			sequence = line[-length:]
			if sequence in sequenceOcurrences:
				sequenceOcurrences[sequence] += 1
			else:
				sequenceOcurrences[sequence] = 0
	sequenceOcurrences = {sequence: occurrences for sequence, occurrences in sequenceOcurrences.items() if occurrences > minimumOfOcurrences}
	print(sequenceOcurrences)

def task3(file):
	
	minimumOfOcurrences = 10000
	hits = 0
	counter = 0
	sequenceOcurrences = dict([])
	# interestingLengths = range(68,78)
	if (file=="s_3_sequence_1M.txt"):
		minimumOfOcurrences = 100000
		interestingLengths = range(40,51)
	elif (file=="s_1-1_1M.txt"):
		minimumOfOcurrences = 300000
		interestingLengths = range(68,78)
	elif (file=="seqset3.txt"):
		minimumOfOcurrences = 50000
		interestingLengths = range(40,52)

	with open(file,'r') as f:
		for line in f: 
			line = line.rstrip()
			for length in interestingLengths:
				sequence = line[-length:]
				if sequence in sequenceOcurrences:
					sequenceOcurrences[sequence] += 1
				else:
					sequenceOcurrences[sequence] = 0
	sequenceOcurrences = {sequence: occurrences for sequence, occurrences in sequenceOcurrences.items() if occurrences > minimumOfOcurrences}
	print(sequenceOcurrences)

	for sequence, _ in sequenceOcurrences.items():
		sequenceInv = sequence[::-1]
		occurrences = task1(file, sequenceInv)
		print(str(sequence) +" has " + str(occurrences)+ " occurrences")

	return hits

def task31(file): #s_1-1_1M.txt
	#file = "s_1-1_1M.txt"
	adapter = "GAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAA"
	adapter= adapter[::-1]
	adapter_tree = SuffixTree()
	adapter_tree.updateTree(adapter)
	remainingLengths = dict([])
	#remainingLength = 76 - len(r)
	hits=0
	with open(file, 'r') as f:
		for line in f:
			line = line.rstrip()
			line = line[::-1]
			stuffInTheTree = adapter_tree.isInTree(line) 
			if (stuffInTheTree != []):
				hits += 1
				suffix = stuffInTheTree[-1]
				remainingLength = 76 - len(suffix)
				if remainingLength not in remainingLengths:
					remainingLengths[remainingLength] = 0
				else:
					remainingLengths[remainingLength] += 1

	del adapter_tree
	print(remainingLengths)
	print(hits) # should output 646364


def readSequence(sequence, adapter, pMatch):
	while(sequence!=""):
		if ( ( sum ( [ sequence[i]==adapter[i] for i in range(len(sequence)) ] ) / len(sequence) ) >= pMatch ):
			break
		else:
			sequence 	= sequence[1:]
			adapter 	= adapter[:-1]

	return sequence

def task22OLD(file, adapter):
	
	# adapter_tree = SuffixTree()
	# adapter_tree.updateTree(adapter)

	hits = 0
	counter = 0
	with open(file, 'r') as f:
		for line in f:

			line = line.rstrip()
			#line = line[::-1]
			# if (adapter_tree.isInTree(line) != []):
			# 	hits += 1
			[tableD, tableL] = functionTables(line, adapter)
			outputArr 		= [] #cambiar a nombre más descriptivo porfiss
			lastColumnD		= tableD[1:,-1]
			lastColumnL		= tableL[1:,-1]
			result			= [[b/(a+b),(a+b)] for a,b in zip(lastColumnD,lastColumnL)]
			result			= result[::-1]
			hit 			= next((index for [value, index] in result if value>=0.9),None)
			if hit:
				hits+=1
			counter=counter+1
			print(counter/1000000)
			#print(hits)
			# guardar este valor para luego
	print(hits) # should output 646364

	
	return hits

def functionTables(text, pattern):

	len_pattern 	= len(pattern)
	len_text 		= len(text)

	tableD = np.zeros((len_pattern+1, len_text+1))
	tableL = np.zeros((len_pattern+1, len_text+1))
	tableD[:,0] = np.arange(len_pattern+1)

	for i,c_pattern in enumerate(pattern,1):
		for j, c_sequence in enumerate(text,1):
			up			= tableD[i-1,j] + 1
			left		= tableD[i,j-1] + 1
			equal		= tableD[i-1, j-1] if (c_sequence == c_pattern) else math.inf
			tableD[i,j]	= min(left, up, equal)

	for i,c_pattern in enumerate(pattern,1):
		for j, c_sequence in enumerate(text,1):
			if (tableD[i,j] == tableD[i-1,j] + 1):
				tableL[i,j] = tableL[i-1,j] 
			elif (tableD[i,j] == tableD[i-1,j-1] + 0 if (c_sequence == c_pattern) else 1):
				tableL[i,j] = tableL[i-1,j-1] + 1
			else:
				tableL[i,j] = tableL[i,j-1] + 1

	# print("tableD:")
	# print(tableD)
	# print("tableL:")
	# print(tableL)

	return [tableD, tableL]

def wordSuffixes(string):
	return [ string[:iterator] for iterator, _ in enumerate(string,1) ]

def buildMatrixD(text, pattern):

	len_pattern 	= len(pattern)
	len_text 		= len(text)

	tableD = np.zeros((len_pattern+1, len_text+1))
	tableL = np.zeros((len_pattern+1, len_text+1))
	tableD[:,0] = np.arange(len_pattern+1)

	for i,c_pattern in enumerate(pattern,1):
		for j, c_sequence in enumerate(text,1):
			up			= tableD[i-1,j] + 1
			left		= tableD[i,j-1] + 1
			equal		= tableD[i-1, j-1] if (c_sequence == c_pattern) else math.inf
			tableD[i,j]	= min(left, up, equal)

	for i,c_pattern in enumerate(pattern,1):
		for j, c_sequence in enumerate(text,1):
			if (tableD[i,j] == tableD[i-1,j] + 1):
				tableL[i,j] = tableL[i-1,j] 
			elif (tableD[i,j] == tableD[i-1,j-1] + 0 if (c_sequence == c_pattern) else 1):
				tableL[i,j] = tableL[i-1,j-1] + 1
			else:
				tableL[i,j] = tableL[i,j-1] + 1

	print("tableD:")
	print(tableD)
	print("tableL:")
	print(tableL)

	return [tableD, tableL]

def buildNextColumn(columnD, columnL, character, pattern):
	
	nextColumnD = np.zeros(len(pattern) + 1)
	nextColumnL = np.zeros(len(pattern) + 1)

	for it, c_pattern in enumerate(pattern,1):

		nextColumnD[it] = nextColumnD[it - 1] + 1
		nextColumnL[it] = nextColumnL[it - 1]
		charactersAreEquals = ( c_pattern ==  character )

		if ( columnD[it -1] + ( not charactersAreEquals ) < nextColumnD[it] ):
			nextColumnD[it] = columnD[it - 1] + (not charactersAreEquals)
			nextColumnL[it] = columnL[it - 1] + 1
		if ( columnD[it] + 1 < nextColumnD[it] ):
			nextColumnD[it] = columnD[it] + 1
			nextColumnL[it] = columnL[it] + 1

	return [nextColumnD, nextColumnL] 

def indexOfMaximumBoundedEntry(orderedColumn, k):
	reverseColumn 	= orderedColumn[::-1]
	auxIndex 		= next(index for value,index in enumerate(reverseColumn) if value>=k)
	return auxIndex



main()
