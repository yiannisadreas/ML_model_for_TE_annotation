import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.model_selection import train_test_split
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from Bio import SeqIO
from torch.utils.data import TensorDataset, DataLoader
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, classification_report, confusion_matrix 
import seaborn as sns
from sklearn.utils import resample



# Set device for GPU/CPU
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")



# Function to one-hot encode a sequence
def one_hot_encode(seq, max_len = 8300):
	mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
	one_hot = np.zeros((4, max_len), dtype=np.float32)

	for i, nucleotide in enumerate(seq[:max_len].upper()):
		if nucleotide in mapping:
			one_hot[mapping[nucleotide], i] = 1.0
	return one_hot




def load_sequences_from_folder(folder_path, label, max_len=8300):
	sequences = []
	labels = []
	for file in os.listdir(folder_path):
		file_path = os.path.join(folder_path, file)
		for record in SeqIO.parse(file_path, "fasta"):
			encoded_seq = one_hot_encode(str(record.seq), max_len)
			sequences.append(encoded_seq)
			labels.append(label)
	return sequences, labels






DNA_folder = r"C:\Users\yiann\Desktop\projects\TE_annotation\DNA_transposons"
LTRs_folder = r"C:\Users\yiann\Desktop\projects\TE_annotation\LTRs"
non_TE_folder = r"C:\Users\yiann\Desktop\projects\TE_annotation\random_sequences"
max_len = 8300
batch_size = 32


LTRs_seq, LTRs_labels = load_sequences_from_folder(LTRs_folder, 2, max_len)
DNA_seq, DNA_labels = load_sequences_from_folder(DNA_folder, 1, max_len)
non_TE_seq, non_TE_labels = load_sequences_from_folder(non_TE_folder, 0, max_len)

min_class_size = min(len(DNA_seq), len(LTRs_seq), len(non_TE_seq))

DNA_seq_sampled, DNA_labels_sampled = resample(
    DNA_seq, DNA_labels,
    replace=False,
    n_samples=min_class_size,
    random_state=42
)
LTRs_seq_sampled, LTRs_labels_sampled = resample(
    LTRs_seq, LTRs_labels,
    replace=False,
    n_samples=min_class_size,
    random_state=42
)
non_TE_seq_sampled, non_TE_labels_sampled = resample(
    non_TE_seq, non_TE_labels,
    replace=False,
    n_samples=min_class_size,
    random_state=42
)

X_data = DNA_seq_sampled + LTRs_seq_sampled + non_TE_seq_sampled
y_data = DNA_labels_sampled + LTRs_labels_sampled + non_TE_labels_sampled

# Split into train / validation / test (70 / 15 / 15).
# The validation split picks the best checkpoint; the test split is only used once at the end.
X_train, X_temp, y_train, y_temp = train_test_split(X_data, y_data, train_size=0.7, test_size=0.3, random_state=42, stratify=y_data)
X_val, X_test, y_val, y_test = train_test_split(X_temp, y_temp, train_size=0.5, test_size=0.5, random_state=42, stratify=y_temp)

# Weighing the classes (inverse frequency in the training split)
train_counts = np.bincount(y_train, minlength=3)
class_weights = torch.tensor(len(y_train) / (3 * train_counts), dtype=torch.float32).to(device)

X_train = torch.tensor(np.array(X_train), dtype=torch.float32)
X_val = torch.tensor(np.array(X_val), dtype=torch.float32)
X_test = torch.tensor(np.array(X_test), dtype=torch.float32)

y_train = torch.tensor(y_train, dtype=torch.long)
y_val = torch.tensor(y_val, dtype=torch.long)
y_test = torch.tensor(y_test, dtype=torch.long)


train_dataset = TensorDataset(X_train, y_train)
val_dataset = TensorDataset(X_val, y_val)
test_dataset = TensorDataset(X_test, y_test)

train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
val_loader = DataLoader(val_dataset, batch_size=batch_size)
test_loader = DataLoader(test_dataset, batch_size=batch_size)

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
	nn.Linear(256, 3)
	# No Softmax here: CrossEntropyLoss expects raw logits. Apply softmax at inference.
).to(device)


loss = nn.CrossEntropyLoss(weight=class_weights)
optimizer = optim.Adam(model.parameters(), lr=0.0001)

best_val_f1 = 0.0
best_model_path = 'best_model.pth'


num_epochs = 500
for epoch in range(num_epochs):
	model.train()
	for batch_X, batch_y in train_loader:
		batch_X = batch_X.to(device)
		batch_y = batch_y.to(device)
		predictions = model(batch_X)
		BCELoss = loss(predictions, batch_y)
		BCELoss.backward()
		optimizer.step()
		optimizer.zero_grad()

	if (epoch + 1) % 100 == 0:
		model.eval()
		with torch.no_grad():
			# Training metrics
			train_preds_list = []
			train_labels_list = []
			for batch_X, batch_y in train_loader:
				batch_X = batch_X.to(device)
				preds = model(batch_X)
				train_preds_list.append(preds.cpu())
				train_labels_list.append(batch_y)
			train_preds = torch.cat(train_preds_list)
			train_labels = torch.cat(train_labels_list)
			train_predicted_labels = torch.argmax(train_preds, dim=1)
			train_accuracy = accuracy_score(train_labels, train_predicted_labels)
			train_precision = precision_score(train_labels, train_predicted_labels, average='macro')
			train_recall = recall_score(train_labels, train_predicted_labels, average='macro')
			train_f1 = f1_score(train_labels, train_predicted_labels, average='macro')

			# Validation metrics
			val_preds_list = []
			val_labels_list = []
			for batch_X, batch_y in val_loader:
				batch_X = batch_X.to(device)
				preds = model(batch_X)
				val_preds_list.append(preds.cpu())
				val_labels_list.append(batch_y)
			val_preds = torch.cat(val_preds_list)
			val_labels = torch.cat(val_labels_list)
			val_predicted_labels = torch.argmax(val_preds, dim=1)
			val_accuracy = accuracy_score(val_labels, val_predicted_labels)
			val_precision = precision_score(val_labels, val_predicted_labels, average='macro')
			val_recall = recall_score(val_labels, val_predicted_labels, average='macro')
			val_f1 = f1_score(val_labels, val_predicted_labels, average='macro')

			# Save the best model based on validation F1 score
			if val_f1 > best_val_f1:
				best_val_f1 = val_f1
				torch.save(model.state_dict(), best_model_path)
				print(f"Saved new best model at epoch {epoch+1} with val_f1: {val_f1:.4f}")

			# Print metrics
			print(
				f'Epoch [{epoch+1}/{num_epochs}], '
				f'Train Acc: {train_accuracy:.4f}, Val Acc: {val_accuracy:.4f}, '
				f'Train Prec: {train_precision:.4f}, Val Prec: {val_precision:.4f}, '
				f'Train Rec: {train_recall:.4f}, Val Rec: {val_recall:.4f}, '
				f'Train F1: {train_f1:.4f}, Val F1: {val_f1:.4f}'
			)

			print(classification_report(val_labels, val_predicted_labels, target_names=["non-TE", "DNA", "LTR"]))

# Final evaluation of the best checkpoint on the held-out test split
model.load_state_dict(torch.load(best_model_path, map_location=device))
model.eval()
test_preds_list = []
test_labels_list = []
with torch.no_grad():
	for batch_X, batch_y in test_loader:
		batch_X = batch_X.to(device)
		preds = model(batch_X)
		test_preds_list.append(preds.cpu())
		test_labels_list.append(batch_y)
test_preds = torch.cat(test_preds_list)
test_labels = torch.cat(test_labels_list)
test_predicted_labels = torch.argmax(test_preds, dim=1)
print(f"Test F1 of best model (val_f1 {best_val_f1:.4f}): {f1_score(test_labels, test_predicted_labels, average='macro'):.4f}")
print(classification_report(test_labels, test_predicted_labels, target_names=["non-TE", "DNA", "LTR"]))

cm = confusion_matrix(test_labels, test_predicted_labels)
sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", xticklabels=["non-TE", "DNA", "LTR"], yticklabels=["non-TE", "DNA", "LTR"])
plt.xlabel("Predicted")
plt.ylabel("True")
plt.show()
