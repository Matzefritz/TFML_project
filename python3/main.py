import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Import Data
path = './../../diterpene_shuf.csv'
data_set = pd.read_csv(path)

row1 = data_set.iloc[0]
frequencies = [row1[f'a{i}'] for i in range(6, 42, 2) ]
multiplicities = [row1[f'a{i}'] for i in range(7, 43, 2) ]

multiplicity_map = {'s': 1, 'd': 2, 't': 3, 'q': 4}

multiplicities_numerical = []
for m in enumerate(multiplicities) :
  multiplicities_numerical.append(multiplicity_map.get(m[1], 0))


#
# --- Plot one diterpene for a start ---
#

plt.style.use('default')

fig, ax = plt.subplots(figsize=(14, 5))

markerline, stemlines, baseline = ax.stem(frequencies, multiplicities_numerical)

# Customize appearance
plt.setp(markerline, marker='o', markersize=6, color='blue')
plt.setp(stemlines, linewidth=2, color='black')

# X-axis ticks and labels
ax.set_xticks(frequencies)
ax.set_xticklabels([f"{f:.1f}" for f in frequencies], rotation=90, fontsize=9)

# Labels and formatting
ax.set_xlabel("Chemical Shift (ppm)", fontsize=12)
ax.set_ylabel("Multiplicity (encoded)", fontsize=12)
ax.set_title("13C NMR Spectrum of First Diterpene", fontsize=14)

ax.invert_xaxis()
ax.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()


#
# --- Heatmap of all the data ---
#

n_samples = len(data_set)
n_peaks = 19
multiplicity_map = {'s': 1, 'd': 2, 't': 3, 'q': 4}

spectra_matrix = np.zeros((n_samples, n_peaks))

for idx, row in data_set.iterrows():
    spectrum = []
    for i in range(7, 44, 2):
        m = row[f'a{i}']
        spectrum.append(multiplicity_map.get(m, 0))
    spectra_matrix[idx, :] = spectrum

# Plot heatmap
plt.figure(figsize=(14, 6))
plt.imshow(spectra_matrix, aspect='auto', cmap='viridis')
plt.colorbar(label="Multiplicity (encoded)")
plt.xlabel("Peak Index (1 to 20)")
plt.ylabel("Spectrum Index")
plt.title("All Diterpene Spectra (Encoded Multiplicities)")
plt.tight_layout()
plt.show()




