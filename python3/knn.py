import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from dataclasses import dataclass, field
from typing import List
from collections import defaultdict

import random
import my_plot

from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split

import import_csv

data = import_csv.import_data()
#
# --- Reconfigure data for Popper ---
#

simple = True

if simple :

	rows, cols = (len(data), 4)
	feature_vector = [[0 for _ in range(cols)] for _ in range(rows)]
	group_vector = [0] * rows

	for i, struct in enumerate(data) :
		count_s = import_csv.count_multiplicity(struct.atoms, 's')
		count_d = import_csv.count_multiplicity(struct.atoms, 'd')
		count_t = import_csv.count_multiplicity(struct.atoms, 't')
		count_q = import_csv.count_multiplicity(struct.atoms, 'q')
		
		feature_vector[i] = [count_s, count_d, count_t, count_q]
		group_vector[i] = struct.group
else :

	multiplicity_map = {'s': 1, 'd': 2, 't': 3, 'q': 4}

	rows, cols = (len(data), 20)
	feature_vector = [[0 for _ in range(cols)] for _ in range(rows)]
	group_vector = [0] * rows

	for i, struct in enumerate(data) :
		for j in range(0,20,2) :
			feature_vector[i][j] = struct.atoms[int(j/2)].frequency
			feature_vector[i][j+1] = multiplicity_map.get(struct.atoms[int(j/2)].multiplicity.lower(), 0)
		group_vector[i] = struct.group


#print(feature_vector)
#print(group_vector)


#
# --- Using SKLEARN for knn ---
#

knn = KNeighborsClassifier(n_neighbors=3)

X_train, X_test, y_train, y_test = train_test_split(feature_vector, group_vector, test_size=0.1)
knn.fit(X_train, y_train)
y_pred = knn.predict(X_test)

print(classification_report(y_test, y_pred))


my_plot.plot_interclass_distance_matrix(feature_vector, group_vector)