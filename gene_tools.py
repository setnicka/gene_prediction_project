#!/usr/bin/python3

start_codon = 'ATG'
stop_codons = ['TAA', 'TAG', 'TGA']


# Function to get all ORF (open reading frames)
# -> format: [
# 	[ startIndex, endIndex, complement = False ]
# ]
#
# WARNING: Internally it uses arrays indexed from 0, but in gene software
# there are used positions indexed from 1 (beware when input/output gene positions)
def GetORF(DNA):
	return(
		_GetORFwithOffset(DNA, 0)
		+ _GetORFwithOffset(DNA, 1)
		+ _GetORFwithOffset(DNA, 2)
		+ _GetORFwithOffset(DNA, 0, True)
		+ _GetORFwithOffset(DNA, 1, True)
		+ _GetORFwithOffset(DNA, 2, True)
	)


# Internal function to find ORF with given direction and offset
def _GetORFwithOffset(DNA, offset, reverse=False):
	if reverse:
		DNA = DnaReverseComplement(DNA)

	started = False
	start = 0
	positions = []
	for i in range(offset, len(DNA) - 2, 3):
		if not started:
			if DNA[i:i + 3] == start_codon:
				started = True
				start = i
		else:
			if DNA[i:i + 3] in stop_codons:
				positions.append([start, i + 2, reverse])
				started = False
	# Append last open segment
	if started:
		positions.append([start, len(DNA) - 1, reverse])

	# In reverse search rotate direction (index 0 = the last one)
	if reverse:
		for item in positions:
			(item[0], item[1]) = (len(DNA) - item[1] - 1, len(DNA) - item[0] - 1)

	return positions


def DnaComplement(dna):
	complement = []
	for c in dna:
		if c == "A":
			complement.append("T")
		elif c == "T":
			complement.append("A")
		elif c == "C":
			complement.append("G")
		elif c == "G":
			complement.append("C")
	return "".join(complement)


def DnaReverse(dna):
	temp = list(dna)
	temp.reverse()
	return "".join(temp)


def DnaReverseComplement(dna):
	return DnaComplement(DnaReverse(dna))


# print(GetORF("CTGCAGACGAAACCTCTTGATGTAGTTGGCCTGACACCGACAATAATGAAGACTACCGTCTTACTAACAC"))
# print(GetORF(DnaReverseComplement("AAATGAA")))
