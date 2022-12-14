import argparse, os, sys, errno, subprocess, csv

'''
TODO - make this script smart => stop when queue is full, move submitted to own directory
'''

def main():
    parser = argparse.ArgumentParser(description="Run submission script.")
    parser.add_argument("--job_dir", type=str,  help="Where are the job submission scripts to be queued?")

    args = parser.parse_args()
    job_dir = args.job_dir

    if not os.path.exists(job_dir):
        print("Unable to find job directory")
        exit(-1)

    job_files = [job_file for job_file in os.listdir(job_dir) if ".sb" in job_file]

    for job_file in job_files:
        job_path = os.path.join(job_dir, job_file)
        print(f"Submitting {job_path}")
        cmd = f"cd {job_dir}; sbatch {job_file}"
        # print(cmd)
        subprocess.run(cmd, shell=True)

if __name__ == "__main__":
    main()