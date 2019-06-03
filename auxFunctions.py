import numpy as np

# Auxiliary functions

def wordSuffixes(text):
		return [ text[:iterator] for iterator, _ in enumerate(text,1) ]

def reverseIndexMaxBound(orderedColumn, k):
	reverseColumn 	= orderedColumn[::-1]
	auxIndex 		= next(index for value,index in enumerate(reverseColumn) if value>=k)
	return auxIndex