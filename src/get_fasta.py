import subprocess
import os

#create bed files from coordinates
def write_bed_file(regions, output_path):
	with open(output_path, "w") as f:
		chrom, start, end = regions
		f.write(f"{chrom}\t{start}\t{end}\n")


#Counts each diferent file
counter = 0

#extract location of TE
TEs = {}
with open("/home/andreas/TE_annotation/D_melanogaster_repeats.out", 'r') as file:
	for line in file:
		parts = line.strip().split()
		if len(parts) < 14:
			continue
		if parts[4].startswith("chr"):
			chromosome = parts[4].replace("chr", "")
			TE_family = parts[10]
			if TE_family.startswith("LTR") or TE_family.startswith("DNA"):
				TE_family_short = TE_family[:3]
				TE_start = parts[5]
				TE_end = parts[6]
				TE = [chromosome, TE_start, TE_end]
				filename = f"{TE_family_short}_{counter:05d}.bed"
				path = os.path.join("/home/andreas/TE_annotation/TEs_bed_files", filename)

				write_bed_file(TE, path)
				counter += 1


genome_path = "/home/andreas/TE_annotation/D_melanogaster_genome.fasta"
fasta_dir = "/home/andreas/TE_annotation/TEs_sequences"

for bed_file in os.listdir("/home/andreas/TE_annotation/TEs_bed_files"):
	bed_file_path = os.path.join("/home/andreas/TE_annotation/TEs_bed_files", bed_file)
	output_file_fasta = os.path.join(fasta_dir, bed_file.replace(".bed", ".fasta"))
	bed_tools_command = [ "bedtools", "getfasta", "-fo", output_file_fasta, "-fi", genome_path, "-bed", bed_file_path]
	result = subprocess.run(bed_tools_command, capture_output=True, text=True)
