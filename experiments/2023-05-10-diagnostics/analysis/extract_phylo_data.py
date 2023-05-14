import sys
import json
import pandas as pd
import glob
from collections import Counter
import time
import csv

dfs = []
dirs = glob.glob(sys.argv[1])


def ancestor_or_descendant(phylo_children, phylo_parents, tax1, tax2):
    curr = tax1
    to_check = []
    parent = phylo_children[curr]
    if parent != -1:
        to_check.append((parent, True))

    if curr in phylo_parents:
        children = phylo_parents[curr]
        for child in children:
            to_check.append((child, False))

    while to_check:
        curr = to_check.pop(0)
        if curr[0] == tax2:
            if curr[1]:
                return 1
            else:
                return -1

        if curr[1]:
            parent = phylo_children[curr[0]]
            if parent != -1:
                to_check.append((parent, True))
        else:
            if curr[0] in phylo_parents:
                children = phylo_parents[curr[0]]
                for child in children:
                    to_check.append((child, False))
    return 0


for dir in dirs:
    df = pd.read_csv(dir+"/run_config.csv", index_col="parameter")
    df = df.T
    print(dir+"/phylo_50000.csv " + str(time.time()), flush=True)
    # phylo = pd.read_csv(dir+"/phylo_50000.csv", memory_map=True, usecols=["id", "ancestor_list", "traits_estimation_source_ids"])
    # phylo.columns = phylo.columns.str.replace(' ', '')
    # phylo["ancestor_list"] = phylo["ancestor_list"].apply(lambda x: int(json.loads(x)[0]) if x.strip().lower() != "[none]" else -1)

    phylo_parents = {}
    phylo_children = {}
    estimates = {}
    with open(dir+"/phylo_50000.csv") as infile:
        header = infile.readline().split(",")
        id_col = header.index("id")
        anc_col = header.index("ancestor_list")
        estimation_sources_col = header.index("traits_estimation_source_ids")

        csv_reader = csv.reader(infile)
        for line in csv_reader:
            id_val = int(line[id_col])
            anc_val = line[anc_col] 
            if anc_val.strip().lower() == "[none]":
                anc_val = -1
            else:
                anc_val = json.loads(anc_val)[0]
            phylo_children[id_val] = anc_val
            estimates[id_val] = json.loads(line[estimation_sources_col])
            
            if anc_val in phylo_parents:
                phylo_parents[anc_val].append(id_val)
            else:
                phylo_parents[anc_val] = [id_val]



    # phylo_parents = phylo.groupby("ancestor_list")["id"].unique().to_dict()

    # phylo.set_index("id", inplace=True)
    # phylo.sort_index(inplace=True)
    # phylo_children = phylo.to_dict()["ancestor_list"]

    ancestor_count = 0
    descendant_count = 0
    self_count = 0
    other_count = 0
    outside_count = 0

    for ind in phylo_children:
        ids = Counter(estimates[ind])
        for id in ids:
            if id == 0:
                continue
            elif id == ind:
                self_count += ids[id]
            elif id not in phylo_children:
                outside_count += ids[id]
            else:
                res = ancestor_or_descendant(phylo_children, phylo_parents, ind, id)
                if res == 1:
                    ancestor_count += ids[id]
                elif res == -1:
                    descendant_count += ids[id]
                else:
                    other_count += ids[id]

    df["ancestor_count"] = ancestor_count
    df["descendant_count"] = descendant_count
    df["self_count"] = self_count
    df["other_count"] = other_count
    df["outside_count"] = outside_count
    dfs.append(df)

all_data = pd.concat(dfs)
all_data.to_csv(sys.argv[2], index=False)