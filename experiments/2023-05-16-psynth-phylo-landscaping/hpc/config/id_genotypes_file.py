import os, argparse

# directory: path to directory to find file(s) in
# contains: list of substrings the candidate files must contain
# identifier_fun: function that takes a list of candidates and returns one or more
def identify_file(directory, contains_substrs, identifier_fun):
    candidates = [cand for cand in os.listdir(directory) if all([req in cand for req in contains_substrs]) ]
    return identifier_fun(candidates)

def id_fun(candidates):
    updates = [int(cand.split("_")[-1].split(".")[0]) for cand in candidates]
    max_update = max(updates)
    return f"phylo_genoetypes_{max_update}.sgp"

def main():
    parser = argparse.ArgumentParser(description="Find genotype file")
    parser.add_argument("--dir", type=str, help="Where is the output directory")
    args = parser.parse_args()
    output_dir = args.dir

    fname = identify_file(
        output_dir,
        [".sgp", "phylo_genotypes_"],
        id_fun
    )

    print(fname)
    exit(0)



if __name__ == "__main__":
    main()