import numpy as np
from auxFunctions import *

class Stack:
	def __init__(self, text="", pattern="", k=0):
		self.text 				= text
		self.pointerText		= 0
		self.pattern	 		= pattern
		self.k					= k
		self.variableM			= len(pattern)
		self.columnLength		= self.variableM +1
		self.wordSuffixes 		= wordSuffixes(text)
		self.availableStates 	= wordSuffixes(text)
		self.eliminatedStates	= [] 	# this might be repetitive tho
		self.initColumnD		= np.zeros(self.columnLength)
		self.initColumnL		= np.zeros(self.columnLength)
		self.columnD			= np.zeros(self.columnLength)
		self.columnL			= np.zeros(self.columnLength)
		self.output				= [] 
		# self.columnsD			
		# self.columnsL


	@property
	def pointedChar(self):
		return self.text[self.pointerText]

	@property 
	def currState(self):
		return self.text[:self.pointerText]
	
	@property
	def previousState(self):
		return self.currState[1:]

	@property
	def getToPreviousState(self):
		self.pointerText = self.pointerText -1
	


def augmentPointer(self):
	self.pointerText += 1

def resetPointer(self):
	self.pointerText = 0


# def initColumn(columnName):
# 	columnName = np.zeros(self.columnLength)
	
def functionTables(text, pattern):

	len_pattern 	= len(pattern)
	len_text 		= len(text)

	tableD = np.zeros((len_pattern+1, len_sequence+1))
	tableL = np.zeros((len_pattern+1, len_sequence+1))
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

	#def dp(tableD, tableL, t):

def buildNextColumns(self):
	
	nextColumnD = np.zeros(self.columnLength)
	nextColumnL = np.zeros(self.columnLength)
	newChar 	= self.pointedChar

	for it, charPattern in enumerate(self.pattern,1):

		nextColumnD[it] 	= nextColumnD[it - 1] + 1
		nextColumnL[it] 	= nextColumnL[it - 1]
		charactersAreEquals = ( charPattern ==  newChar )

		if ( self.columnD[it -1] + ( not charactersAreEquals ) < nextColumnD[it] ):
			nextColumnD[it] = self.columnD[it - 1] + (not charactersAreEquals)
			nextColumnL[it] = self.columnL[it - 1] + 1
		
		if ( self.columnD[it] + 1 < nextColumnD[it] ):
			nextColumnD[it] = self.columnD[it] + 1
			nextColumnL[it] = self.columnL[it] + 1

	self.columnD = nextColumnD
	self.columnL = nextColumnL


def algorithmC(text, pattern, k):
	states = '0' + [ text[:iterator] for iterator, _ in enumerate(text,1) ]
	eliminatedStates = ['0']
	columnD = np.zeros(len(pattern) + 1)
	columnL = np.zeros(len(pattern) + 1)
	search(0, columnD, columnL)


def search(self, r, columnD, columnL):

	# define state = new_char+oldtext
	if self.currState in self.availableStates:
		[nextColumnD, nextColumnL] = buildNextColumns(self.columnD, self.columnL, self.pointedChar)
		indexOfMaxBoundEntry 	= self.columnLength - indexOfMaximumBoundedEntry(self.columnD, self.k)
		viablePrefixLength 		= self.columnL[indexOfMaximumBoundedEntry]

		if len(self.currState)> viablePrefixLength:
			while(len(self.currState)>viablePrefixLength and state in states):
				eliminatedStates.append(state)
				self.currState = self.previousState
		if len(state)==viablePrefixLength and state in states:
			if d[indexM]<=k:
				saveToStack(state)
			eliminatedStates.append(state)
			auxStateW = state
			precedentOfW = auxStateW[1:]
			# while(precedentOfW and ):
			# search(state, self.newChar)
			# append w' to eliminatedStates for all w' such that 

def saveToStack(self,state):
	self.output.append(state)

# def search(root, columnD, columnL):
# 	for each state s=g(r,a) Do
# 		if (not eliminated(s)):
# 			[nextColumnD, nextColumnL] = buildNextColumn(columnD, columnL, a) 



def indexOfMaximumBoundedEntry(orderedColumn, k):
	reverseColumn 	= orderedColumn[::-1]
	auxIndex 		= next(index for value,index in enumerate(reverseColumn) if value>=k)
	return auxIndex

	# IDEA FOR EX3 AND 4: SUFFIX TREE WITH DIFFERENT LAYERS. 