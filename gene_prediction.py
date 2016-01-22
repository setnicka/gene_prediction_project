#!/usr/bin/python3

import math
import operator
from gene_settings import *
import gene_tools

import cg_islands
import statistical_prediction


# Some statistics:
def ComputeStats(sequences):
	print("Sekvence:\t\t\t\t# genů\tPrůměrná délka genu\t# ORF\tPrůměrná délka ORF")
	for seq in sequences:
		real_genes = gene_tools.ReadGenes(sequences_dir + seq + ".gb.genes")
		fasta = gene_tools.ReadFasta(sequences_dir + seq + ".fasta")
		DNA = list(fasta.values())[0]
		ORF = gene_tools.GetORF(DNA)
		real_sum = 0
		for gene in real_genes:
			real_sum += gene[1] - gene[0] + 1
		ORF_sum = 0
		for gene in ORF:
			ORF_sum += gene[1] - gene[0] + 1

		print("%s\t%s\t%s\t%s\t%s" % (seq.ljust(35), len(real_genes), real_sum / len(real_genes), len(ORF), ORF_sum / len(ORF)))


def Evaluate(real_genes, found_genes):
	real_genes.sort(key=operator.itemgetter(0, 1))
	found_genes.sort(key=operator.itemgetter(0, 1))

	found = 0
	notfound = 0
	misfound = 0

	index_real = 0
	index_found = 0
	while index_real < len(real_genes) and index_found < len(found_genes):
		r = real_genes[index_real]
		f = found_genes[index_found]
		if r[0] < f[0] or (r[0] == f[0] and r[0] < f[0]):
			index_real += 1
			notfound += 1
		elif r[0] > f[0] or (r[0] == f[0] and r[0] > f[0]):
			index_found += 1
			misfound += 1
		else:
			# Same start and ends
			if r[0] != f[0]:
				"WRONG"
			found += 1
			index_found += 1
			index_real += 1

	return (found, notfound, misfound)


# Function for testing different gene prediction methods and pretty printing
def Test(method, params, sequences):
	print("Sekvence:\t\t\t\tFound\t!Found\tMisfound")
	for seq in sequences:
		real_genes = gene_tools.ReadGenes(sequences_dir + seq + ".gb.genes")
		fasta = gene_tools.ReadFasta(sequences_dir + seq + ".fasta")
		DNA = list(fasta.values())[0]
		ORF = gene_tools.GetORF(DNA)

		found_genes = method(DNA, ORF, *params)

		print(seq.ljust(35), end="\t")
		print("%d\t%d\t%d" % Evaluate(real_genes, found_genes))


################################################################################

def LengthPrediction(fasta, ORF, treshold):
	return [gene for gene in ORF if gene[1] - gene[0] >= treshold]


print("=====\nLength prediction:")
Test(LengthPrediction, (150, ), sequences)

################################################################################

statistic_treshold = 0.295 ** (statistical_prediction.model_pre + statistical_prediction.model_post)
length_treshold = 150

print("=====\nStatistical prediction:")
model = statistical_prediction.StatisticalPredictionTrain(training_sequences)
Test(statistical_prediction.StatisticalPrediction, (model, statistic_treshold, length_treshold), sequences)

################################################################################

print("=====\nCG islands prediction:")
log_matrix = cg_islands.CGislandsPredictionTrain(training_sequences)
Test(cg_islands.CGislandsPrediction, (log_matrix, 0.004, length_treshold), sequences)
