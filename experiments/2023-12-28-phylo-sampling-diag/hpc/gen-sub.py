'''
Generate slurm job submission scripts - one per condition
'''

import argparse, os, sys, pathlib
# from asyncio import base_events
# from email.policy import default
from pyvarco import CombinationCollector
sys.path.append(os.path.join(pathlib.Path(os.path.dirname(os.path.abspath(__file__))).parents[2], "scripts"))
import utilities as utils

default_seed_offset = 1000
default_account = "devolab"
default_num_replicates = 30
default_job_time_request = "120:00:00"
default_job_mem_request = "4G"
default_total_generations = 300000

job_name = "12-28"
executable = "diagnostics"

base_script_filename = "./base_script.txt"

# Create combo object to collect all conditions we'll run
combos = CombinationCollector()

fixed_parameters = {
    "POP_SIZE": "500",
    "STOP_MODE": "generations",
    "MAX_GENS": "50000",
    "MAX_EVALS": "500000000",
    "INIT_POP_RAND": "1",
    "GENE_LOWER_BND": "0",
    "GENE_UPPER_BND": "1",
    "ACCURACY": "0.98",
    "TARGET": "100",
    "CREDIT": "0",
    "DIAGNOSTIC_DIMENSIONALITY": "100",
    "MUTATE_PER_SITE_RATE": "0.007",
    "MUTATE_MEAN": "0",
    "MUTATE_STD": "1",
    "TOURNAMENT_SIZE": "7",
    "OUTPUT_SUMMARY_DATA_INTERVAL": "100",
    "PRINT_INTERVAL": "20",
    "SNAPSHOT_INTERVAL": "10000",
    "EVAL_MAX_PHYLO_SEARCH_DEPTH": "8"
}

special_decorators = ["__DYNAMIC", "__COPY_OVER"]

combos.register_var("eval__COPY_OVER")
combos.register_var("DIAGNOSTIC")
combos.register_var("SELECTION")

combos.add_val(
    "eval__COPY_OVER",
    [

        "-EVAL_MODE down-sample -TEST_DOWNSAMPLE_RATE 0.1 -EVAL_FIT_EST_MODE none",
        "-EVAL_MODE down-sample -TEST_DOWNSAMPLE_RATE 0.1 -EVAL_FIT_EST_MODE ancestor",
        "-EVAL_MODE indiv-rand-sample -TEST_DOWNSAMPLE_RATE 0.1 -EVAL_FIT_EST_MODE ancestor",
        "-EVAL_MODE phylo-informed-sample -TEST_DOWNSAMPLE_RATE 0.1 -EVAL_FIT_EST_MODE ancestor",

        "-EVAL_MODE down-sample -TEST_DOWNSAMPLE_RATE 0.01 -EVAL_FIT_EST_MODE none",
        "-EVAL_MODE down-sample -TEST_DOWNSAMPLE_RATE 0.01 -EVAL_FIT_EST_MODE ancestor",
        "-EVAL_MODE indiv-rand-sample -TEST_DOWNSAMPLE_RATE 0.01 -EVAL_FIT_EST_MODE ancestor",
        "-EVAL_MODE phylo-informed-sample -TEST_DOWNSAMPLE_RATE 0.01 -EVAL_FIT_EST_MODE ancestor",

        "-EVAL_MODE full -TEST_DOWNSAMPLE_RATE 1.0 -NUM_COHORTS 0 -EVAL_FIT_EST_MODE none",
    ]
)

combos.add_val(
    "SELECTION",
    [
        "lexicase",
        "tournament"
    ]
)

combos.add_val(
    "DIAGNOSTIC",
    [
        "exploitation-rate",
        "contradictory-objectives",
        "multipath-exploration"
    ]
)

def main():
    parser = argparse.ArgumentParser(description="Run submission script.")
    parser.add_argument("--data_dir", type=str, help="Where is the output directory for phase one of each run?")
    parser.add_argument("--config_dir", type=str, help="Where is the configuration directory for experiment?")
    parser.add_argument("--repo_dir", type=str, help="Where is the repository for this experiment?")
    parser.add_argument("--job_dir", type=str, default=None, help="Where to output these job files? If none, put in 'jobs' directory inside of the data_dir")
    parser.add_argument("--replicates", type=int, default=default_num_replicates, help="How many replicates should we run of each condition?")
    parser.add_argument("--seed_offset", type=int, default=default_seed_offset, help="Value to offset random number seeds by")
    parser.add_argument("--account", type=str, default=default_account, help="Value to use for the slurm ACCOUNT")
    parser.add_argument("--time_request", type=str, default=default_job_time_request, help="How long to request for each job on hpc?")
    parser.add_argument("--mem", type=str, default=default_job_mem_request, help="How much memory to request for each job?")

    # Load in command line arguments
    args = parser.parse_args()
    data_dir = args.data_dir
    config_dir = args.config_dir
    job_dir = args.job_dir
    repo_dir = args.repo_dir
    num_replicates = args.replicates
    hpc_account = args.account
    seed_offset = args.seed_offset
    job_time_request = args.time_request
    job_memory_request = args.mem
    # patch_mode = args.patch

    # Load in the base slurm file
    base_sub_script = ""
    with open(base_script_filename, 'r') as fp:
        base_sub_script = fp.read()

    # Get list of all combinations to run
    combo_list = combos.get_combos()

    # Calculate how many jobs we have, and what the last id will be
    num_jobs = num_replicates * len(combo_list)
    print(f'Generating {num_jobs} across {len(combo_list)} files!')
    print(f' - Data directory: {data_dir}')
    print(f' - Config directory: {config_dir}')
    print(f' - Repository directory: {repo_dir}')
    print(f' - Replicates: {num_replicates}')
    print(f' - Account: {hpc_account}')
    print(f' - Time Request: {job_time_request}')
    print(f' - Seed offset: {seed_offset}')

    # If no job_dir provided, default to data_dir/jobs
    if job_dir == None:
        job_dir = os.path.join(data_dir, "jobs")

    # Create job file for each condition
    cur_job_id = 0
    cond_i = 0
    generated_files = set()
    for condition_dict in combo_list:
        cur_seed = seed_offset + (cur_job_id * num_replicates)
        filename_prefix = f'RUN_C{cond_i}'
        file_str = base_sub_script
        file_str = file_str.replace("<<TIME_REQUEST>>", job_time_request)
        file_str = file_str.replace("<<MEMORY_REQUEST>>", job_memory_request)
        file_str = file_str.replace("<<JOB_NAME>>", job_name)
        file_str = file_str.replace("<<CONFIG_DIR>>", config_dir)
        file_str = file_str.replace("<<REPO_DIR>>", repo_dir)
        file_str = file_str.replace("<<EXEC>>", executable)
        file_str = file_str.replace("<<JOB_SEED_OFFSET>>", str(cur_seed))
        file_str = file_str.replace("<<ACCOUNT_NAME>>", hpc_account)

        ###################################################################
        # Configure the run
        ###################################################################
        file_str = file_str.replace("<<RUN_DIR>>", \
            os.path.join(data_dir, f'{filename_prefix}_'+'${SEED}'))

        # Format commandline arguments for the run
        run_param_info = {key:condition_dict[key] for key in condition_dict if not any([dec in key for dec in special_decorators])}
        # Add fixed paramters
        for param in fixed_parameters:
            if param in run_param_info: continue
            run_param_info[param] = fixed_parameters[param]
        # Set random number seed
        run_param_info["SEED"] = '${SEED}'

        ###################################################################
        # Build avida commandline parameters string
        ###################################################################
        fields = list(run_param_info.keys())
        fields.sort()
        set_params = [f"-{field} {run_param_info[field]}" for field in fields]
        copy_params = [condition_dict[key] for key in condition_dict if "__COPY_OVER" in key]
        run_params = " ".join(set_params + copy_params)
        ###################################################################

        # Add run commands to run the experiment
        cfg_run_commands = ''
        # Set the run
        cfg_run_commands += f'RUN_PARAMS="{run_params}"\n'

        # By default, add all commands to submission file.
        array_id_run_info = {
            array_id: {
                "experiment": True
            }
            for array_id in range(1, num_replicates+1)
        }
        array_id_to_seed = {array_id:(cur_seed + (array_id - 1)) for array_id in array_id_run_info}

        # Track which array ids need to be included. If none, don't need to output this file.
        active_array_ids = []
        inactive_array_ids = []
        run_sub_logic = ""
        # NOTE - this is setup to (fairly) easily incorporate job patching logic, but not currently incorporated
        for array_id in range(1, num_replicates+1):
            # If this run is totally done, make note and continue.
            # if not any([array_id_run_info[array_id][field] for field in array_id_run_info[array_id]]):
            #     inactive_array_ids.append(array_id)
            #     continue
            # This run is not done already. Make note.
            active_array_ids.append(array_id)

            run_logic = "if [[ ${SLURM_ARRAY_TASK_ID} -eq "+str(array_id)+" ]] ; then\n"

            # (1) Run experiment executable
            run_commands = ''
            run_commands += 'echo "./${EXEC} ${RUN_PARAMS}" > cmd.log\n'
            run_commands += './${EXEC} ${RUN_PARAMS} > run.log\n'

            run_logic += run_commands
            # run_logic += analysis_commands
            run_logic += "fi\n\n"
            run_sub_logic += run_logic

        # -- Set the SLURM array id range parameter --
        array_id_range_param = ""
        if len(active_array_ids) == num_replicates:
            array_id_range_param = f"1-{num_replicates}"
        else:
            array_id_range_param = ",".join([str(array_id) for array_id in active_array_ids])

        # -- add run commands to file str --
        file_str = file_str.replace("<<ARRAY_ID_RANGE>>", array_id_range_param)
        file_str = file_str.replace("<<CFG_RUN_COMMANDS>>", cfg_run_commands)
        file_str = file_str.replace("<<RUN_COMMANDS>>", run_sub_logic)

        ###################################################################
        # Write job submission file (if any of the array ids are active)
        ###################################################################
        # Report active/inactive
        # print(f"RUN_C{cond_i}:")
        # print(f" - Active: " + ", ".join([f"RUN_C{cond_i}_{array_id_to_seed[array_id]}" for array_id in active_array_ids]))
        # print(f" - Inactive: " + ", ".join([f"RUN_C{cond_i}_{array_id_to_seed[array_id]}" for array_id in inactive_array_ids]))

        if len(active_array_ids):
            utils.mkdir_p(job_dir)
            with open(os.path.join(job_dir, f'{filename_prefix}.sb'), 'w') as fp:
                fp.write(file_str)

        # Update condition id and current job id
        cur_job_id += 1
        cond_i += 1

if __name__ == "__main__":
    main()
