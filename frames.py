#!/usr/bin/python3

# Functions to get all ORF (open reading frames)
# -> format: [
# 	[ startIndex, endIndex, complement = False ]
# ]
#
# WARNING: Internally it uses arrays indexed from 0, but in gene software
# there are used positions indexed from 1 (beware when input/output gene positions)

start_codon = 'ATG'
stop_codons = ['TAA', 'TAG', 'TGA']


def GetORF(DNA):
	return(
		GetORFwithOffset(DNA, 0)
		+ GetORFwithOffset(DNA, 1)
		+ GetORFwithOffset(DNA, 2)
		+ GetORFwithOffset(DNA, 0, True)
		+ GetORFwithOffset(DNA, 1, True)
		+ GetORFwithOffset(DNA, 2, True)
	)


def GetORFwithOffset(DNA, offset, reverse=False):
	if reverse:
		DNA = DnaReverseComplement(DNA)

	started = False
	start = 0
	positions = []
	for i in range(offset, len(DNA) - 2, 3):
		print(i, DNA[i:i + 3])
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
		positions.append([start, len(DNA), reverse])
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


print(GetORF("CTGCAGACGAAACCTCTTGATGTAGTTGGCCTGACACCGACAATAATGAAGACTACCGTCTTACTAACAC"))
