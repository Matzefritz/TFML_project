import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from dataclasses import dataclass, field
from typing import List
from collections import defaultdict
import random

from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split	


@dataclass
class Atom:
	frequency: float
	multiplicity: str

@dataclass
class Structure:
	name: str
	group: str
	atoms: List[Atom] = field(default_factory=list)


def count_multiplicity(atoms, string):
	count = 0
	for atom in atoms :
		if atom.multiplicity == string :
			count = count + 1

	return count


#
# --- Import Data
#

def import_data() :

	path = './../../diterpene_shuf.csv'
	data_set = pd.read_csv(path)

	data: List[Structure] = []
	for i, row in data_set.iterrows():
		row = data_set.iloc[i]
		atom_list: List[Atom] = []
		for j in range(5,44,2):
			mult = row[f'a{j}']
			freq = row[f'a{j+1}']
			atom_list.append(Atom(freq, mult))

		data.append(Structure(row['a4'], row['a45c'], atom_list))



	row1 = data_set.iloc[0]
	frequencies = [row1[f'a{i}'] for i in range(6, 42, 2) ]
	multiplicities = [row1[f'a{i}'] for i in range(7, 43, 2) ]

	multiplicity_map = {'s': 1, 'd': 2, 't': 3, 'q': 4}

	multiplicities_numerical = []
	for m in enumerate(multiplicities) :
	  multiplicities_numerical.append(multiplicity_map.get(m[1], 0))


	#
	# --- Add rules for reduction from the paper ---
	#

	data_reduced = data.copy()

	reducing = False

	if reducing :
		for i, struct in enumerate(data):
			for j, atom in enumerate(data[i].atoms):
				m = data[i].atoms[j].multiplicity
				f = data[i].atoms[j].frequency

				if m == 's':
					if 64.5 <= f <= 95:
						data_reduced[i].atoms[j].multiplicity = 'd'
					elif 96 <= f <= 114:
						data_reduced[i].atoms[j].multiplicity = 't'
					elif 115 <= f <= 165:
						data_reduced[i].atoms[j].multiplicity = 'd'
					elif 165 <= f <= 188:
						data_reduced[i].atoms[j].multiplicity = 'q'
					elif f > 188:
						data_reduced[i].atoms[j].multiplicity = 't'

				elif m == 'd':
					if 64.5 <= f <= 95:
						data_reduced[i].atoms[j].multiplicity = 't'
					elif 96 <= f <= 104:
						data_reduced[i].atoms[j].multiplicity = 'q'
					elif 105 <= f <= 180:
						data_reduced[i].atoms[j].multiplicity = 't'
					elif f > 180:
						data_reduced[i].atoms[j].multiplicity = 'q'

				elif m == 't':
					if 59 <= f <= 90:
						data_reduced[i].atoms[j].multiplicity = 'q'
					elif f > 90:
						data_reduced[i].atoms[j].multiplicity = 'q'

		data = data_reduced

	return data