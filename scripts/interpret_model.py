import matplotlib.pyplot as plt
from captum.attr import IntegratedGradients
import torch
import numpy as np
import os
from Bio import SeqIO
import torch.nn as nn


device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

class_names = {0: "non-TE", 1: "DNA", 2: "LTR"}

# Ensure compatibility with KMP # Short-term fix and might cause silent bugs
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

# Function to one-hot encode a sequence
def one_hot_encode(seq, max_len=8300):
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    one_hot = np.zeros((4, max_len), dtype=np.float32)

    for i, nucleotide in enumerate(seq[:max_len].upper()):
        if i >= max_len:
            break
        if nucleotide in mapping:
            one_hot[mapping[nucleotide], i] = 1.0
    return one_hot


   
    

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
	nn.Linear(256, 3),
	nn.Softmax(dim=1)
).to(device)


model.load_state_dict(torch.load("best_model.pth", map_location=device))
model.eval()

def load_fasta_sequences(file_path):
    records = []
    for record in SeqIO.parse(file_path, "fasta"):
        records.append((record.id, str(record.seq).upper()))
    return records 

# Compute Integrated Gradients
def compute_ig(model, input_tensor, target_class):
    input_tensor = input_tensor.unsqueeze(0).to(device)
    input_tensor.requires_grad_()
    ig = IntegratedGradients(model)
    attributions, delta = ig.attribute(input_tensor, target=target_class, return_convergence_delta=True)
    return attributions.squeeze(0).detach().cpu().numpy()

# Figure plotting function
def plot_ig_figure(attributions, sequence, k=100):
    importance = attributions.sum(axis=0)
    sequence = np.array(list(sequence))
    
    top_k_indices = np.argsort(importance)[-k:][::-1] # Descending order
    
    top_k_importances = importance[top_k_indices]
    top_k_bases = sequence[top_k_indices]
    top_k_positions = top_k_indices
    
    labels = [f"{base} ({pos})" for base, pos in zip(top_k_bases, top_k_positions)]
    
    print("Min:", np.min(importance))
    print("Max:", np.max(importance))
    print("Mean:", np.mean(importance))
    print("Sum:", np.sum(importance))
    print(top_k_positions)
    
    plt.figure(figsize=(12, 4))
    plt.bar(labels, top_k_importances, color='blue', alpha=0.7)
    plt.xticks(rotation=45, ha='right', fontsize=9)
    plt.xlabel("Position")
    plt.ylabel("Importance")
    plt.title(f"Top {k} Most Important Bases (Integrated Gradients)")
    plt.tight_layout()
    plt.show()
    
if __name__ == "__main__":
    fasta_file = r"C:\Users\yiann\Desktop\projects\TE_annotation\DNA_transposons\DNA_00259.fasta"  
    records = load_fasta_sequences(fasta_file)
    seq_id, dna_sequence = records[0]
    
    encoded_seq = one_hot_encode(dna_sequence)
    input_tensor = torch.tensor(encoded_seq, dtype=torch.float32)

    with torch.no_grad():
        output = model(input_tensor.unsqueeze(0).to(device))
        predicted_class = torch.argmax(output, dim=1).item()
        confidence = output[0, predicted_class].item()
        print(f"Predicted class: {class_names[predicted_class]} (confidence: {confidence:.4f})")


    attributions = compute_ig(model, input_tensor, target_class=predicted_class)
    plot_ig_figure(attributions, dna_sequence)