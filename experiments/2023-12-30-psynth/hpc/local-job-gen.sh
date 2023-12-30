#!/usr/bin/env bash

REPLICATES=50
EXP_SLUG=2023-12-30-psynth
ACCOUNT=devolab
SEED_OFFSET=60000
JOB_TIME=48:00:00
JOB_MEM=16G
PROJECT_NAME=phylogeny-informed-subsampling

SCRATCH_EXP_DIR=./test/data/${PROJECT_NAME}
REPO_DIR=/Users/lalejina/devo_ws/${PROJECT_NAME}
HOME_EXP_DIR=${REPO_DIR}/experiments

DATA_DIR=${SCRATCH_EXP_DIR}/${EXP_SLUG}
JOB_DIR=${SCRATCH_EXP_DIR}/${EXP_SLUG}/jobs
CONFIG_DIR=${HOME_EXP_DIR}/${EXP_SLUG}/hpc/config

python3 gen-sub.py --runs_per_subdir 500 --time_request ${JOB_TIME} --mem ${JOB_MEM} --data_dir ${DATA_DIR} --config_dir ${CONFIG_DIR} --repo_dir ${REPO_DIR} --replicates ${REPLICATES} --job_dir ${JOB_DIR} --account ${ACCOUNT} --seed_offset ${SEED_OFFSET}