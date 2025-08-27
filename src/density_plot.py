import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import sys


arguments = sys.argv
region_path = arguments[1]
#region_file = "/home/andreas/TEs_Project/RepeatMasker/wntA/Bicyclus_anynana.bed"
for region_file in os.listdir(region_path):
	if not region_file.endswith("closest_gene.bed"):
		continue
	with open(os.path.join(region_path, region_file), 'r') as file:
		lines = file.readlines()
		up_region_parts = lines[0].strip().split()
		chromosome = up_region_parts[0]
		print(up_region_parts)
		up_gene_start = int(up_region_parts[6])
		up_gene_end = int(up_region_parts[7])
		down_region_parts = lines[1].strip().split()
		down_gene_start = int(down_region_parts[6])
		down_gene_end = int(down_region_parts[7])
#		region_start = int((down_gene_start + down_gene_end) / 2)
#		region_end = int((up_gene_start + up_gene_end) / 2)
#		print(region_file)


	name = region_file.replace(".bed", "")
	name = name.replace("closest_gene", "mapping")
	gene_name = arguments[2]

	gene_path = arguments[3] #C:\Users\yiann\Desktop\projects\TE_annotation\_50_ QueryCov\wntA genomes
	for gene_file in os.listdir(gene_path):
		if not gene_file.endswith(".bed"):
			continue
		with open(os.path.join(gene_path, gene_file), 'r') as file:
			line = file.read()
			gene_parts = line.strip().split()
			if chromosome == gene_parts[0]:
				gene_region = (int(gene_parts[1]), int(gene_parts[2]))
				print(gene_region)
				print(gene_file)
				species = gene_file.replace(".bed", "")
				break



	repeats_path = r"C:\Users\yiann\Desktop\projects\TE_annotation\TEs_predictions"
#	repeats_file = "/home/andreas/TEs_Project/Repeats_IGV/wntA/Bicyclus_anynana.bed"
	insertions = {}
	for repeats_file in os.listdir(repeats_path):
		if not repeats_file.endswith(".bed"):
			continue
		if species not in repeats_file:
			continue
		with open(os.path.join(repeats_path, repeats_file), 'r') as file:
			for line in file:
				repeat_parts = line.strip().split()
				if chromosome != repeat_parts[0]:
					continue
				repeat_start = int(repeat_parts[1])
				if repeat_start > up_gene_start or repeat_start < down_gene_end:
					continue
				repeat_end = int(repeat_parts[2])
				TE_type = repeat_parts[3]
				repeat_cor = [repeat_start, repeat_end]
				if TE_type not in insertions:
					insertions[TE_type] = []
				insertions[TE_type].append(repeat_cor)
				repeats_file_name = repeats_file




	# GFF file for the gene
	gff_file_path = r"C:\Users\yiann\Desktop\projects\TE_annotation\Genome_annotations"
	locus_mapping = {}
	for gff_file in os.listdir(gff_file_path):
		if not gff_file.endswith(".gff"):
			continue
		if species not in gff_file:
			continue
		species_gff = os.path.join(gff_file_path, gff_file)
		with open(species_gff, 'r') as gff_file:
			for line in gff_file:
				gff_parts = line.strip().split()
				if gff_parts[0] != chromosome:
					continue
				if gff_parts[2] == "exon" or gff_parts[2] == "tRNA" or gff_parts[2] == "lnc_RNA":
					element = gff_parts[2]
					if element not in locus_mapping:
						locus_mapping[element] = []
					element_start = int(gff_parts[3])
					element_end = int(gff_parts[4])
					if element_start >= down_gene_start and element_end <= up_gene_end:
						element_region = (element_start, element_end)
						locus_mapping[element].append(element_region)
			


	fig, ax = plt.subplots(figsize=(24, 9))
	ax.set_xlim(down_gene_start, up_gene_end)
	ax.set_ylim(0.5, 4.5)
	ax.set_xlabel("Genomic Position (bp)")
	ax.set_yticks([0.8, 1.8, 2.8, 3.8])
	ax.set_yticklabels(["Gene: " + gene_name, "Elements", "DNA Transposons", "LTR Retrotransposons"])
	ax.set_title("Genome Locus: Gene and Nearby TE Insertions")

	ax.add_patch(patches.Rectangle((gene_region[0], 1 - 0.4), gene_region[1] - gene_region[0], 0.4, facecolor="navy", edgecolor="black"))
	ax.add_patch(patches.Rectangle((down_gene_start, 1 - 0.4), down_gene_end - down_gene_start, 0.4, facecolor="red", edgecolor="black"))
	ax.add_patch(patches.Rectangle((up_gene_end, 1 - 0.4), up_gene_start - up_gene_end, 0.4, facecolor="red", edgecolor="black"))
	
 
	for element, regions in locus_mapping.items():
		for start, end in regions:
			if element == "exon":
				ax.add_patch(patches.Rectangle((start, 2 - 0.4), end - start, 0.4, facecolor="lightgreen", edgecolor="black", label="Elements"))
			elif element == "tRNA":
				ax.add_patch(patches.Rectangle((start, 2 - 0.4), end - start, 0.4, facecolor="yellow", edgecolor="black", label="Elements"))
			elif element == "lnc_RNA":
				ax.add_patch(patches.Rectangle((start, 2 - 0.4), end - start, 0.4, facecolor="darkmagenta", edgecolor="black", label="Elements"))
			
	
 
	if "DNA" in insertions:
		for start, end in insertions["DNA"]:
			ax.add_patch(patches.Rectangle((start, 3 - 0.4), end - start, 0.4, facecolor= "orange", edgecolor="black", label="DNA Transposon"))

	if "LTR" in insertions:
		for start, end in insertions["LTR"]:
			ax.add_patch(patches.Rectangle((start, 4 - 0.4), end - start, 0.4, facecolor= "cyan", edgecolor="black", label="LTR Retrotransposon"))

	plt.tight_layout()
	plt.savefig(name + ".png")
	plt.close(fig)
#	plt.show()