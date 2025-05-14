import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from collections import defaultdict
import seaborn as sns

def plot_interclass_distance_matrix(X, y, metric='euclidean'):
    """
    X: list of feature vectors
    y: list of class labels
    metric: distance metric (default 'euclidean')
    """
    # Group feature vectors by class
    class_to_vecs = defaultdict(list)
    for xi, label in zip(X, y):
        class_to_vecs[label].append(xi)

    classes = sorted(class_to_vecs.keys())
    dist_matrix = np.zeros((len(classes), len(classes)))

    # Compute mean pairwise distances between all class pairs
    for i, class_i in enumerate(classes):
        for j, class_j in enumerate(classes):
            vecs_i = np.array(class_to_vecs[class_i])
            vecs_j = np.array(class_to_vecs[class_j])

            # Pairwise distances between two groups
            dists = cdist(vecs_i, vecs_j, metric=metric)
            dist_matrix[i, j] = np.mean(dists)

    # Plot heatmap
    plt.figure(figsize=(12, 10))
    ax = sns.heatmap(dist_matrix, annot=True, fmt=".2f", cmap="viridis",
                     xticklabels=classes, yticklabels=classes, square=True, cbar=True)

    plt.title(f"Average Pairwise Distances Between Classes ({metric})")
    plt.xlabel("Class")
    plt.ylabel("Class")
    plt.tight_layout()
    plt.show()
