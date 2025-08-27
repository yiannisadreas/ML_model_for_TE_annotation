import torch
import torch.nn as nn
import numpy as np
from Bio import SeqIO
import os


# Define device
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


class_names = {0: "non-TE", 1: "DNA", 2: "LTR"}

# Recreate the model architecture (must match training)
model = nn.Sequential(
    nn.Conv1d(in_channels=4, out_channels=32, kernel_size=16),
    nn.ReLU(),
    nn.MaxPool1d(kernel_size=2),
    nn.Conv1d(in_channels=32, out_channels=64, kernel_size=8),
    nn.ReLU(),
    nn.MaxPool1d(kernel_size=2),
    nn.Conv1d(in_channels=64, out_channels=128, kernel_size=8),
    nn.ReLU(),
    nn.MaxPool1d(kernel_size=2),
    nn.Conv1d(in_channels=128, out_channels=256, kernel_size=6),
    nn.ReLU(),
    nn.MaxPool1d(kernel_size=16),
    nn.AdaptiveMaxPool1d(1),
    nn.Flatten(),
    nn.Dropout(0.3),
    nn.Linear(256,3),
    nn.Softmax(dim=1)
).to(device)

# Load model
model.load_state_dict(torch.load("best_model.pth", map_location=device))
model.eval()


def load_fasta_sequences(file_path):
    records = []
    for record in SeqIO.parse(file_path, "fasta"):
        records.append((record.id, str(record.seq).upper()))
    return records 



def sliding_windows(seq, window_size=8300, step_size=500):
    for i in range(0, len(seq) - window_size + 1, step_size):
        yield i, seq[i:i+window_size]



def one_hot_encode(seq, max_len = 8300):
	mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
	one_hot = np.zeros((4, max_len), dtype=np.float32)

	for i, nucleotide in enumerate(seq[:max_len].upper()):
		if nucleotide in mapping:
			one_hot[mapping[nucleotide], i] = 1.0
	return one_hot


for file in os.listdir(r"C:\Users\yiann\Desktop\projects\TE_annotation\lepidoptera_genomes"):
    if not file.endswith(".fasta"):
        continue
    file_path = os.path.join(r"C:\Users\yiann\Desktop\projects\TE_annotation\lepidoptera_genomes", file)
    print(f"Processing {file_path}...")

    # Load the model
    predicted_regions = []

    for seq_name, full_seq in load_fasta_sequences(file_path):
        for start, window_seq in sliding_windows(full_seq, window_size=8300):
            input_tensor = torch.tensor(one_hot_encode(window_seq)).unsqueeze(0).to(device)
            
            with torch.no_grad():
                output = model(input_tensor)
                predicted_class = torch.argmax(output, dim=1).item()
                confidence = output[0, predicted_class].item()
                if predicted_class in [1, 2] and confidence > 0.8:  # Threshold for confidence
                    predicted_regions.append((seq_name, start, start + 8300, class_names[predicted_class], confidence))

    filename = file.replace(".fasta", "")
    tsv_file_path = os.path.join(r"C:\Users\yiann\Desktop\projects\TE_annotation\TEs_predictions", f"{filename}_TEs_predictions.tsv")

    # Write the merged regions to a file
    with open(tsv_file_path, "w") as f:
        f.write("seq_name\tstart\tend\tclass\tconfidence\n")
        for region in predicted_regions:
            f.write(f"{region[0]}\t{region[1]}\t{region[2]}\t{region[3]}\t{region[4]:.2f}\n")




"""
def merge_regions(regions):
    if not regions:
        return []
    regions.sort()  # Sort by seq_name and start
    merged = [regions[0]]
    for curr in regions[1:]:
        last = merged[-1]
        # If same seq_name and overlapping/adjacent
        if curr[0] == last[0] and curr[1] <= last[2]:
            # Merge by extending the end
            merged[-1] = (last[0], last[1], max(last[2], curr[2]))
        else:
            merged.append(curr)
    return merged
"""




