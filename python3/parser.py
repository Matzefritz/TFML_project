import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from dataclasses import dataclass, field
from typing import List
from collections import defaultdict
import random

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

print(data[0].atoms)

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
reducing = True

if reducing :

	data_reduced = data.copy()

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


#
# --- Reconfigure data for Popper ---
#


# Background File
with open("./../popper/my_parsed_pl_files/popper_run/bk.pl", "w") as f :
	for i, struct in enumerate(data):
		f.write('prop(' + struct.name + "," + str(count_multiplicity(struct.atoms, 's')) + "," + str(count_multiplicity(struct.atoms, 'd')) + "," + str(count_multiplicity(struct.atoms, 't')) + "," + str(count_multiplicity(struct.atoms, 'q')) + ").\n")


# Examples File

# Groud data after group names
grouped_data = defaultdict(list)
for struct in data :
	grouped_data[struct.group].append(struct)


for group, structs in grouped_data.items() :
	if len(structs) >= 10 :
		# Sample 5 structs from the selected group
		true_samples = random.sample(structs, 10)

		# Sample 5
		flag = True
		while flag :
			false_samples = random.sample(data, 10)
			flag = False
			for sam in false_samples :
				if sam.group == group :
					flag = True

	# Write an examples file for each group
	with open(f'./../popper/my_parsed_pl_files/examples_{group}.pl', "w") as f :
		for struct in true_samples :
			f.write('pos(belongs(' + struct.name + ",g" + struct.group + ")).\n")
		for struct in false_samples :
			f.write('neg(belongs(' + struct.name + ",g" + true_samples[0].group + ")).\n")


# Bias file
#	with open("./../popper/my_parsed_pl_files/popper_run/bias.pl", "w") as f :
#		f.write("head_pred(belongs,2).\n")
#		f.write("body_pred(prop,5).\n")
#		f.write("type(prop,(mol_id,int,int,int,int)).")




"""
%–– what we’re learning
head_pred(belongs,2).

%–– background predicates
body_pred(prop,5).
body_pred(gt,2).
body_pred(lt,2).

%–– types
type(belongs,(molecule,atom)).
type(prop,(molecule,int,int,int,int)).
type(gt,(num,num)).
type(lt,(num,num)).

%–– which B‐preds may appear in a belongs/2 clause
determination(belongs/2,prop/5).
determination(belongs/2,gt/2).
determination(belongs/2,lt/2).

%–– ensure exactly one “molecule” variable in each clause
:-
    clause(C),
    #count{V : var_type(C,V,molecule)} != 1.

%–– FOIL’s “one rule, ≤6 literals” bias
max_clauses(1).
max_body(6).
	



"""