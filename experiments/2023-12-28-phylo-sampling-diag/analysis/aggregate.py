'''
Aggregate data
'''

import argparse, os, sys, pathlib
sys.path.append(
    os.path.join(
        pathlib.Path(
            os.path.dirname(os.path.abspath(__file__))).parents[2],
            "scripts"
    )
)
import utilities as utils

run_identifier = "RUN_"

run_cfg_fields = [
    "SEED",
    "POP_SIZE",
    "STOP_MODE",
    "MAX_GENS",
    "MAX_EVALS",
    "INIT_POP_RAND",
    "GENE_LOWER_BND",
    "GENE_UPPER_BND",
    "ACCURACY",
    "TARGET",
    "CREDIT",
    "DIAGNOSTIC_DIMENSIONALITY",
    "MUTATE_PER_SITE_RATE",
    "MUTATE_MEAN",
    "MUTATE_STD",
    "TOURNAMENT_SIZE",
    "EVAL_MODE",
    "TEST_DOWNSAMPLE_RATE",
    "EVAL_FIT_EST_MODE",
    "NUM_COHORTS",
    "DIAGNOSTIC",
    "SELECTION",
    "EVAL_MAX_PHYLO_SEARCH_DEPTH"
]

agg_exclude_fields = [
    "genome",
    "true_phenotype",
    "evaluated_tests",
    "evaluated_phenotype"
]

time_series_summary_fields = [
    "update",
    "evaluations",
    "pop_optimal_trait_coverage",
    "max_agg_score",
    "num_unique_selected",
    "entropy_selected_ids",
    "parents_optimal_trait_coverage",
    "optimal_trait_coverage_loss",
    "trait_estimations"
]
time_series_phylo_fields = [
    "mean_genotype_pairwise_distance",
    "mrca_changes",
    "mean_genoetype_deleterious_steps",
    "genotype_current_phylogenetic_diversity"

]
time_series_sys_fields = [
    "mrca_depth",
    "diversity",
    "ave_depth"
]
time_series_elite_fields = [
    "eval_true_agg_score_diff",
    "lineage_length",
    "mut_count",
    "deleterious_steps"
]

time_series_cfg_fields = [
    "SEED",
    "DIAGNOSTIC",
    "SELECTION",
    "EVAL_MODE",
    "TEST_DOWNSAMPLE_RATE",
    "EVAL_FIT_EST_MODE",
    "NUM_COHORTS"
]

def main():
    parser = argparse.ArgumentParser(description="Run submission script.")
    parser.add_argument("--data_dir", type=str, help="Where is the base output directory for each run?")
    parser.add_argument("--dump", type=str, help="Where to dump this?", default=".")
    parser.add_argument("--units", type=str, default="interval", choices=["interval", "total"], help="Unit for resolution of time series")
    parser.add_argument("--resolution", type=int, default=1, help="What resolution should we collect time series data at?")

    # Parse command line arguments
    args = parser.parse_args()
    data_dir = args.data_dir
    dump_dir = args.dump
    time_series_units = args.units
    time_series_resolution = args.resolution

    # Verify that the given data directory exits
    if not os.path.exists(data_dir):
        print("Unable to find data directory.")
        exit(-1)

    # Verify time series resolution >= 1
    if time_series_resolution < 1:
        print("Time series resolution must be >= 1")
        exit(-1)

    # Create the directory to dump aggregated data (if it doesn't already exist)
    utils.mkdir_p(dump_dir)

    # Aggregate run directories.
    run_dirs = [run_dir for run_dir in os.listdir(data_dir) if run_identifier in run_dir]
    print(f"Found {len(run_dirs)} run directories.")

    # Create file to hold time series data
    time_series_content = []    # This will hold all the lines to write out for a single run; written out for each run.
    time_series_header = None   # Holds the time series file header (verified for consistency across runs)
    time_series_fpath = os.path.join(dump_dir, f"time_series.csv")

    with open(time_series_fpath, "w") as fp:
        fp.write("")

    # Create data structure to hold summary information for each run (1 element per run)
    summary_header = None
    summary_content_lines = []

    incomplete_runs = []

    # Loop over runs, aggregating data from each.
    total_runs = len(run_dirs)
    cur_run_i = 0
    for run_dir in run_dirs:
        run_path = os.path.join(data_dir, run_dir)
        summary_info = {}                   # Hold summary information about run. Indexed by field.
        time_series_info = {}               # Hold time series information. Indexed by update.

        cur_run_i += 1
        print(f"Processing ({cur_run_i}/{total_runs}): {run_path}")

        # Skip over (but make note of) incomplete runs.
        required_files = [
            os.path.join("output", "summary.csv"),
            os.path.join("output", "elite.csv"),
            os.path.join("output", "phylodiversity.csv"),
            os.path.join("output", "systematics.csv"),
            os.path.join("output", "run_config.csv"),
            os.path.join("cmd.log")
        ]
        incomplete = any([not os.path.exists(os.path.join(run_path, req)) for req in required_files])
        if incomplete:
            print("  - Failed to find all required files!")
            incomplete_runs.append(run_dir)
            continue

        ############################################################
        # Extract commandline configuration settings (from cmd.log file)
        # cmd_log_path = os.path.join(run_path, "cmd.log")
        # cmd_params = utils.extract_params_cmd_log(cmd_log_path, "diagnostics")
        # print(cmd_params)
        # for field in cmd_params:
        #     summary_info[field] = cmd_params[field]

        # Extract configs from run_config.csv file
        cfg_path = os.path.join(run_path, "output", "run_config.csv")
        cfg_data = utils.read_csv(cfg_path)
        cmd_params = {}
        for line in cfg_data:
            param = line["parameter"]
            value = line["value"]
            cmd_params[param] = value
            if param in run_cfg_fields:
                summary_info[param] = value
        print("Run configuration:", summary_info)
        ############################################################

        ############################################################
        # Extract summary data
        run_data_path = os.path.join(run_path, "output", "summary.csv")
        run_data = utils.read_csv(run_data_path)

        # Identify final generation
        generations = [int(line["update"]) for line in run_data]
        max_gen = max(generations)
        # Identify evaluations
        evaluations = [int(line["evaluations"]) for line in run_data]
        max_evals = max(evaluations)

        # Isolate data from final generation
        final_data = [line for line in run_data if int(line["update"]) == max_gen]
        assert len(final_data) == 1
        final_data = final_data[0]

        # Add final data to run summary info
        for field in final_data:
            if field in agg_exclude_fields or field in summary_info:
                continue
            summary_info[f"{field}"] = final_data[field]

        # Process time series
        filtered_ts_data = utils.filter_ordered_data(
            run_data,
            time_series_units,
            time_series_resolution
        )
        ts_gens_included = [int(line["update"]) for line in filtered_ts_data]
        ts_gens_included.sort()
        gen_to_ts_step = {ts_gens_included[i]:i for i in range(0, len(ts_gens_included))}
        time_series_info = {gen:{} for gen in ts_gens_included}

        for line in filtered_ts_data:
            gen = int(line["update"])
            time_series_info[gen]["ts_step"] = gen_to_ts_step[gen]
            for field in time_series_summary_fields:
                time_series_info[gen][field] = line[field]

        # Add config info to time series entries
        for gen in time_series_info:
            for field in time_series_cfg_fields:
                time_series_info[gen][field] = summary_info[field]
        ############################################################


        ############################################################
        # Add some logic to detect if this run is incomplete based on config
        stop_condition = cmd_params["STOP_MODE"]

        completed = False
        if stop_condition == "generations":
            completed = int(cmd_params["MAX_GENS"]) <= max_gen
        elif stop_condition == "evaluations":
            completed = int(cmd_params["MAX_EVALS"]) <= max_evals
        else:
            print("Unrecognized stop mode")

        if not completed:
            print("  - Run failed to finish.")
            incomplete_runs.append(run_dir)
            # continue
        ############################################################

        ############################################################
        # Extract elite data
        ############################################################
        elite_data_path = os.path.join(run_path, "output", "elite.csv")
        elite_data = utils.read_csv(elite_data_path)
        # Identify final generation
        generations = [int(line["update"]) for line in elite_data]
        max_gen = max(generations)
        # Identify evaluations
        evaluations = [int(line["evaluations"]) for line in elite_data]
        max_evals = max(evaluations)

        # Isolate data from final generation
        final_data = [line for line in elite_data if int(line["update"]) == max_gen]
        assert len(final_data) == 1
        final_data = final_data[0]

        # Add final data to run summary info
        for field in final_data:
            if field in agg_exclude_fields or field in summary_info:
                continue
            summary_info[f"elite_{field}"] = final_data[field]

        # Process time series
        if completed:
            filtered_ts_data = utils.filter_ordered_data(
                elite_data,
                time_series_units,
                time_series_resolution
            )

            for line in filtered_ts_data:
                gen = int(line["update"])
                time_series_info[gen]["ts_step"] = gen_to_ts_step[gen]
                for field in time_series_elite_fields:
                    time_series_info[gen][field] = line[field]
        ############################################################


        ############################################################
        # Extract phylo  data
        ############################################################
        phylo_data_path = os.path.join(run_path, "output", "phylodiversity.csv")
        phylo_data = utils.read_csv(phylo_data_path)
        # Identify final generation
        generations = [int(line["update"]) for line in phylo_data]
        max_gen = max(generations)
        # Identify evaluations
        evaluations = [int(line["evaluations"]) for line in phylo_data]
        max_evals = max(evaluations)

        # Isolate data from final generation
        final_data = [line for line in phylo_data if int(line["update"]) == max_gen]
        assert len(final_data) == 1
        final_data = final_data[0]

        # Add final data to run summary info
        for field in final_data:
            if field in agg_exclude_fields or field in summary_info:
                continue
            summary_info[f"phylo_{field}"] = final_data[field]

        # Process time series
        if completed:
            filtered_ts_data = utils.filter_ordered_data(
                phylo_data,
                time_series_units,
                time_series_resolution
            )

            for line in filtered_ts_data:
                gen = int(line["update"])
                time_series_info[gen]["ts_step"] = gen_to_ts_step[gen]
                for field in time_series_phylo_fields:
                    time_series_info[gen][field] = line[field]

        ############################################################
        # Extract sys data
        ############################################################
        sys_data_path = os.path.join(run_path, "output", "systematics.csv")
        sys_data = utils.read_csv(sys_data_path)
        # Identify final generation
        generations = [int(line["update"]) for line in sys_data]
        max_gen = max(generations)

        # Isolate data from final generation
        final_data = [line for line in sys_data if int(line["update"]) == max_gen]
        assert len(final_data) == 1
        final_data = final_data[0]

        # Add final data to run summary info
        for field in final_data:
            if field in agg_exclude_fields or field in summary_info:
                continue
            summary_info[f"sys_{field}"] = final_data[field]

        # Process time series
        if completed:
            filtered_ts_data = utils.filter_ordered_data(
                sys_data,
                time_series_units,
                time_series_resolution
            )

            for line in filtered_ts_data:
                gen = int(line["update"])
                time_series_info[gen]["ts_step"] = gen_to_ts_step[gen]
                for field in time_series_sys_fields:
                    time_series_info[gen][field] = line[field]
        ############################################################

        ############################################################
        # Output time series data for this run
        # Compute time series header from time_series_info
        if completed:
            time_series_fields = list(time_series_info[ts_gens_included[0]].keys())
            time_series_fields.sort()
            # If we haven't written the header, write it.
            write_header = False
            if time_series_header == None:
                write_header = True
                time_series_header = ",".join(time_series_fields)
            elif time_series_header != ",".join(time_series_fields):
                print("Time series header mismatch!")
                exit(-1)

            # Write time series content line-by-line
            time_series_content = []
            for u in ts_gens_included:
                time_series_content.append(",".join([str(time_series_info[u][field]) for field in time_series_fields]))
            with open(time_series_fpath, "a") as fp:
                if write_header: fp.write(time_series_header)
                fp.write("\n")
                fp.write("\n".join(time_series_content))
            time_series_content = []
        ############################################################

        ############################################################
        # Add summary_info to aggregate content
        summary_fields = list(summary_info.keys())
        summary_fields.sort()
        if summary_header == None:
            summary_header = summary_fields
        elif summary_header != summary_fields:
            print("Header mismatch!")
            exit(-1)
        summary_line = [str(summary_info[field]) for field in summary_fields]
        summary_content_lines.append(",".join(summary_line))
        ############################################################

    # write out aggregate data
    with open(os.path.join(dump_dir, "aggregate.csv"), "w") as fp:
        out_content = ",".join(summary_header) + "\n" + "\n".join(summary_content_lines)
        fp.write(out_content)

    # Write out incomplete runs, sort them!
    incomplete_runs.sort()
    with open(os.path.join(dump_dir, "incomplete_runs_agg.log"), "w") as fp:
        fp.write("\n".join(incomplete_runs))

if __name__ == "__main__":
    main()


