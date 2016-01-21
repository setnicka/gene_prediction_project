#!/usr/bin/python3

start_codons = ('ATG', 'GTG', 'TTG')
stop_codons = ('TAA', 'TAG', 'TGA')


# Return a dictionary of sequences from FASTA file with the name <filename>.
# The keys are the sequence identifiers, the values are the sequences
# stored as strings.
def ReadFasta(filename, verbose=False):
	dictionary = {}
	identifier = ""
	dna = ""

	def Save():
		if dna and identifier:
			dictionary[identifier] = dna
			if verbose:
				print("The sequence %s with %d characters was read." % (identifier, len(dna)))

	inputFile = open(filename, 'rU')
	try:
		for line in inputFile:
			if line[0] == '>':
				Save()
				identifier = line.split()[0]
				dna = ""
			else:
				dna += line.replace('\n', '')
		Save()
	finally:
		inputFile.close()
	return dictionary


# Parse genes file and return sequence of genes:
# format: [
# 	[ startIndex, endIndex, complement = False ]
# ]
#
# WARNING: Internally it uses arrays indexed from 0, but in gene software
# there are used positions indexed from 1 (beware when input/output gene positions)
# (so when reading sub 1 from each index)
def ReadGenes(filename):
	genes = []

	inputFile = open(filename, 'rU')
	try:
		for line in inputFile:
			line = line.strip()
			complement = False
			if line.startswith("complement"):
				complement = True
				line = line[11:-1]  # Remove the word "complement(" and ending ")"
			start, end = line.split("..")
			start = int(start.replace(">", "").replace("<", ""))
			end = int(end.replace(">", "").replace("<", ""))
			genes.append([start - 1, end - 1, complement])
	finally:
		inputFile.close()
	return genes


# Function to get all ORF (open reading frames)
# format: [
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

	starts = []
	positions = []
	for i in range(offset, len(DNA) - 2, 3):
		if DNA[i:i + 3] in start_codons:
			starts.append(i)
			start = i

		if len(starts) and DNA[i:i + 3] in stop_codons:
			for start in starts:
				positions.append([start, i + 2, reverse])
				starts = []
	# Append last open segment
	if len(starts):
		for start in starts:
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
