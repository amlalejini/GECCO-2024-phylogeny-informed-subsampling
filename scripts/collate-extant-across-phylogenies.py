#!/usr/bin/python3

import glob
import multiprocessing
import os
import re

import pandas as pd


def collate_run(run: str) -> pd.DataFrame:
  print(f"collating run {run}")

  phylogeny_dfs = []
  #for phylogeny_path in glob.glob(f"{run}/output/phylo_*.csv"):
  for phylogeny_path in [
    f"{run}/output/phylo_1000.csv",
    f"{run}/output/phylo_20000.csv",
  ]:
    print(f"reading {phylogeny_path}")
    phylogeny_df = pd.read_csv(phylogeny_path)
    phylogeny_filename = os.path.basename(phylogeny_path)
    update = int(re.search(r'phylo_(\d+)\.csv', phylogeny_filename).group(1))

    phylogeny_df["update"] = update
    phylogeny_dfs.append(
      phylogeny_df[phylogeny_df["destruction_time"] == float("inf")].copy()
    )


  config_df = pd.read_csv(f"{run}/output/run_config.csv")
  config_dict = dict(zip(config_df["parameter"], config_df["value"]))

  run_df = pd.concat(phylogeny_dfs).reset_index(drop=True)

  run_df["EVAL_MODE"] = config_dict["EVAL_MODE"]
  run_df["TEST_DOWN_SAMPLE_RATE"] = config_dict["TEST_DOWNSAMPLE_RATE"]
  run_df["DIAGNOSTIC"] = config_dict["DIAGNOSTIC"]
  run_df["EVAL_FIT_EST_MODE"] = config_dict["EVAL_FIT_EST_MODE"]

  run_df["run"] = run

  return run_df


if __name__ == "__main__":
  runs = glob.glob("RUN_*")
  # note: in practice, requires bugfix Python3.8+
  # https://stackoverflow.com/a/53368158
  with multiprocessing.Pool(16) as pool:
    run_dfs = pool.map(collate_run, runs)

  res_df = pd.concat(run_dfs).reset_index(drop=True)
  res_df.to_csv(
    "~/2023-03-18-diagnostics-phylo-collated.csv.gz",
    index=False,
  )

