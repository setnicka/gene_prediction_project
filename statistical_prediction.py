#!/usr/bin/python3

import gene_tools
from gene_settings import sequences_dir

model_pre = 15  # To catch Pribnow and Gilbert box
model_post = 25  # To catch ribosomal binding site


def StatisticalPredictionTrain(sequences):
	model = None
	for seq in sequences:
		real_genes = gene_tools.ReadGenes(sequences_dir + seq + ".gb.genes")
		fasta = gene_tools.ReadFasta(sequences_dir + seq + ".fasta")
		DNA = list(fasta.values())[0]
		model = _StatisticalPredictionTrain(DNA, real_genes, model)
	return model


def _StatisticalPredictionTrain(fasta, real_genes, model=None):
	if model is None:
		model = {}
		for i in range(-model_pre, model_post):
			model[i] = {'A': 0, 'C': 0, 'T': 0, 'G': 0, 'total': 0}

	complement_fasta = gene_tools.DnaComplement(fasta)

	for gene in real_genes:
		for i in range(-model_pre, model_post):
			if gene[2]:
				index = gene[1] - i
			else:
				index = gene[0] + i

			if index > len(fasta) or index < 0:
				continue
			else:
				if gene[2]:
					base = complement_fasta[index]
				else:
					base = fasta[index]
				model[i][base] += 1
				model[i]['total'] += 1
	return model


def StatisticalPrediction(fasta, ORF, model, statistic_treshold, length_treshold):
	ORF = [gene for gene in ORF if gene[1] - gene[0] >= length_treshold]

	complement_fasta = gene_tools.DnaComplement(fasta)

	output = []
	# For each gene compute score
	for gene in ORF:
		score = 1.0
		for i in range(-model_pre, model_post):
			if gene[2]:
				index = gene[1] - i
			else:
				index = gene[0] + i

			if index > len(fasta) or index < 0:
				continue
			else:
				if gene[2]:
					base = complement_fasta[index]
				else:
					base = fasta[index]

				if base in ('A', 'C', 'T', 'G'):
					score *= (model[i][base] / model[i]['total'])
				else:
					score *= (min(model[i].values()) / model[i]['total'])
		# print(gene, score)
		if score > statistic_treshold:
			output.append(gene)
	return output
