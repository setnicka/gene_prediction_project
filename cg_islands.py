#!/usr/bin/python3

import operator
import math
import gene_tools
from gene_settings import sequences_dir


def CGislandsPredictionTrain(sequences):
	positive_matrix = {}
	negative_matrix = {}
	for base in 'ACTG':
		positive_matrix[base] = {}
		negative_matrix[base] = {}
		for base2 in 'ACTG':
			positive_matrix[base][base2] = 0
			negative_matrix[base][base2] = 0

	for seq in sequences:
		real_genes = gene_tools.ReadGenes(sequences_dir + seq + ".gb.genes")
		fasta = gene_tools.ReadFasta(sequences_dir + seq + ".fasta")
		DNA = list(fasta.values())[0]

		_CGislandsPredictionTrain(DNA, real_genes, positive_matrix, negative_matrix)

	log_matrix = {}
	for base in 'ACTG':
		positive_sum = sum(positive_matrix[base].values())
		negative_sum = sum(negative_matrix[base].values())
		log_matrix[base] = {}
		for base2 in 'ACTG':
			positive_matrix[base][base2] /= positive_sum
			negative_matrix[base][base2] /= negative_sum
			log_matrix[base][base2] = math.log(positive_matrix[base][base2] / negative_matrix[base][base2])
	return log_matrix


def _CGislandsPredictionTrain(fasta, real_genes, positive_matrix, negative_matrix):
	actions = []
	for gene in real_genes:
		actions.append([gene[0], 1])
		actions.append([gene[1], -1])
	actions.sort(key=operator.itemgetter(0))

	in_gene = 0
	actions_index = 0
	for i in range(len(fasta) - 1):
		if actions[actions_index][0] == i:
			in_gene += actions[actions_index][1]
			actions_index += 1
		if in_gene:
			positive_matrix[fasta[i]][fasta[i + 1]] += 1
		else:
			negative_matrix[fasta[i]][fasta[i + 1]] += 1


def CGislandsPrediction(fasta, ORF, log_matrix, CG_treshold, length_treshold):
	ORF = [gene for gene in ORF if gene[1] - gene[0] >= length_treshold]

	output = []
	# For each gene compute CG probability
	for gene in ORF:
		score = 0
		for i in range(gene[0], gene[1]):
			if fasta[i] in 'ACTG' and fasta[i + 1] in 'ACTG':
				score += log_matrix[fasta[i]][fasta[i + 1]]

		if score > CG_treshold * (gene[1] - gene[0]):
			output.append(gene)
	return output
