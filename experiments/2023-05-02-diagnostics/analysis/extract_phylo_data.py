import ALifeStdDev.phylogeny as phylodev
import sys
import ast
import pandas as pd
import glob
from copy import deepcopy
from collections import Counter

dfs = []
dirs = glob.glob(sys.argv[1])

for dir in dirs:
    phylo_files = glob.glob(dir+"/phylo_*.csv")
    df = pd.read_csv(dir+"/run_config.csv", index_col="parameter")
    df = df.T
    # df.reset_index(inplace=True)
    # df.drop("index", inplace=True)
    for f in phylo_files:
        print(f)
        phylo_df = deepcopy(df)
        phylo = phylodev.load_phylogeny_to_networkx(f)
        ancestor_count = 0
        descendant_count = 0
        self_count = 0
        other_count = 0
        outside_count = 0
        for node in phylo.nodes():
            ids = Counter(ast.literal_eval(phylo.nodes[node]['traits_estimation_source_ids']))
            for id in ids:
                if id == 0:
                    continue
                elif id == node:
                    self_count += ids[id]
                elif id not in phylo:
                    outside_count += ids[id]
                elif phylodev.is_ancestor_asexual(phylo, node, id):
                    ancestor_count += ids[id]
                elif phylodev.is_ancestor_asexual(phylo, id, node):
                    descendant_count += ids[id]
                else:
                    other_count += ids[id]
        phylo_df["ancestor_count"] = ancestor_count
        phylo_df["descendant_count"] = descendant_count
        phylo_df["self_count"] = self_count
        phylo_df["other_count"] = other_count
        phylo_df["outside_count"] = outside_count
        phylo_df["time"] = f.split("/")[-1].strip("phylo_.csv")
        dfs.append(phylo_df)

all_data = pd.concat(dfs)
all_data.to_csv(sys.argv[2], index=False)