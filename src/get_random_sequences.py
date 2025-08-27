import random
import os
import subprocess


#create bed files from coordinates
def write_bed_file(regions, output_path):
	with open(output_path, "w") as f:
		chrom, start, end = regions
		f.write(f"{chrom}\t{start}\t{end}\n")


#Dictionary for the TE coordinates inside our genome
TE_coordinates = {}
chromosomes = []
for bed_file in os.listdir("/home/andreas/TE_annotation/TEs_bed_files"):
	bed_file_path = os.path.join("/home/andreas/TE_annotation/TEs_bed_files", bed_file)
	with open(bed_file_path, 'r') as file:
		line = file.readline()
		parts = line.strip().split()
		chromosome = parts[0]
		if len(chromosome) > 5:
			continue
		start = parts[1]
		end = parts[2]
		if chromosome not in TE_coordinates:
			TE_coordinates[chromosome] = []
			chromosomes.append(chromosome)
		else:
			TE_coordinates[chromosome].append((start, end))


chromosomes.remove("L")


#~40000000 the maximum of bases per chromosome
#~3700 files of transposable elements

genome_path = "/home/andreas/TE_annotation/D_melanogaster_genome.fasta"
fasta_dir = "/home/andreas/TE_annotation/random_sequences"


k = 0
i = 0
while i <= 3700:
	random_chromosome = chromosomes[random.randint(0, (len(chromosomes)-1))]
	random_start = random.randint(0, 40000000)
	random_end = random.randint(0, 40000000)
	start = min(random_start, random_end)
	end = max(random_start, random_end)

	if abs(start - end) > 1000 and abs(start - end) < 8000:
		non_overlapping = True
		for existing_start, existing_end in TE_coordinates.get(random_chromosome, []):
			existing_start = int(existing_start)
			existing_end = int(existing_end)
			if not (end <= existing_start or start >= existing_end):
				non_overlapping = False
				break

		if non_overlapping:
			random_range = [random_chromosome, start, end]
			bed_path = os.path.join("/home/andreas/TE_annotation/random_bed_files", f"sequence({i}).bed")
			write_bed_file(random_range, bed_path)
			output_file_fasta = os.path.join(fasta_dir, f"sequence({i}).fasta")
			bed_tools_command = [ "bedtools", "getfasta", "-fo", output_file_fasta, "-fi", genome_path, "-bed", bed_path]
			result = subprocess.run(bed_tools_command, capture_output=True, text=True)
			if result.returncode != 0:
				print("bedtools failed:", result.stderr)
#				continue
			if os.path.getsize(output_file_fasta) == 0:
				os.remove(output_file_fasta)
				continue

			i += 1



















