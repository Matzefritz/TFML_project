import pandas as pd
import random
from dataclasses import dataclass, field
from typing import List
from collections import defaultdict
from copy import deepcopy

# --- your Atom / Structure definitions ---
@dataclass
class Atom:
    frequency: float
    multiplicity: str

@dataclass
class Structure:
    id: str
    group: str
    atoms: List[Atom] = field(default_factory=list)

def count_multiplicity(atoms: List[Atom], m: str) -> int:
    return sum(1 for atom in atoms if atom.multiplicity == m)

# --- 1. Load & parse CSV ---
path = './../../diterpene_shuf.csv'
df = pd.read_csv(path)

data: List[Structure] = []
for _, row in df.iterrows():
    atoms = [ Atom(row[f'a{j}'], row[f'a{j+1}']) for j in range(6, 42, 2) ]
    data.append(Structure(id=row['a4'], group=row['a45c'], atoms=atoms))

# --- 2. Deep-copy & apply reduction rules (from the paper) ---
data_reduced = deepcopy(data)
for struct in data_reduced:
    for atom in struct.atoms:
        m, f = atom.multiplicity, atom.frequency
        new_m = m
        if m == 's':
            if 64.5 <= f <= 95:   new_m = 'd'
            elif 96  <= f <= 114: new_m = 't'
            elif 115 <= f <= 165: new_m = 'd'
            elif 165 <= f <= 188: new_m = 'q'
            elif f > 188:         new_m = 't'
        elif m == 'd':
            if 64.5 <= f <= 95:   new_m = 't'
            elif 96  <= f <= 104: new_m = 'q'
            elif 105 <= f <= 180: new_m = 't'
            elif f > 180:         new_m = 'q'
        elif m == 't':
            if f >= 59:           new_m = 'q'
        # else leave new_m = m
        atom.multiplicity = new_m

# --- 3. Write the single shared background file (bk.pl) ---
bk_path = './../popper/my_parsed_pl_files/popper_run/bk2.pl'
with open(bk_path, 'w') as f:
    for struct in data_reduced:
        s = count_multiplicity(struct.atoms, 's')
        d = count_multiplicity(struct.atoms, 'd')
        t = count_multiplicity(struct.atoms, 't')
        q = count_multiplicity(struct.atoms, 'q')
        f.write(f"prop({struct.id},{s},{d},{t},{q}).\n")
    # Add numeric comparisons:
    f.write("\n")
    f.write("gt(X,Y) :- X > Y.\n")
    f.write("geq(X,Y) :- X >= Y.\n")
    f.write("lt(X,Y) :- X < Y.\n")
    f.write("leq(X,Y) :- X =< Y.\n")

# --- 4. Generate one examples_<group>.pl per class with ≥10 members ---
out_dir = './../popper/my_parsed_pl_files/'
grouped = defaultdict(list)
for struct in data_reduced:
    grouped[struct.group].append(struct)

for group, members in grouped.items():
    if len(members) < 10:
        continue

    pos_samples = random.sample(members, 10)
    # negatives: sample from all *other* groups
    neg_pool = [s for s in data_reduced if s.group != group]
    neg_samples = random.sample(neg_pool, 10)

    ex_path = f"{out_dir}examples_{group}.pl"
    with open(ex_path, 'w') as f:
        for s in pos_samples:
            f.write(f"pos(belongs({s.id},'{group}')).\n")
        for s in neg_samples:
            f.write(f"neg(belongs({s.id},'{group}')).\n")
