#!/usr/bin/env bash

REPLICATES=20
EXP_SLUG=2023-04-02-psynth-small-ds
ACCOUNT=devolab
SEED_OFFSET=6000
JOB_TIME=96:00:00
JOB_MEM=16G
PROJECT_NAME=phylogeny-informed-evaluation

SCRATCH_EXP_DIR=/mnt/scratch/lalejini/data/${PROJECT_NAME}
REPO_DIR=/mnt/home/lalejini/devo_ws/${PROJECT_NAME}
HOME_EXP_DIR=${REPO_DIR}/experiments

DATA_DIR=${SCRATCH_EXP_DIR}/${EXP_SLUG}
JOB_DIR=${SCRATCH_EXP_DIR}/${EXP_SLUG}/jobs
CONFIG_DIR=${HOME_EXP_DIR}/${EXP_SLUG}/hpc/config

python3 gen-sub.py --time_request ${JOB_TIME} --mem ${JOB_MEM} --data_dir ${DATA_DIR} --config_dir ${CONFIG_DIR} --repo_dir ${REPO_DIR} --replicates ${REPLICATES} --job_dir ${JOB_DIR} --account ${ACCOUNT} --seed_offset ${SEED_OFFSET}