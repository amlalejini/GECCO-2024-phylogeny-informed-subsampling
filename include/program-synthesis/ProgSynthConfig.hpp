#pragma once

#include "emp/config/config.hpp"

namespace psynth {

EMP_BUILD_CONFIG(ProgSynthConfig,
  GROUP(WORLD, "How should the world be setup?"),
  VALUE(SEED, int, 0, "Random number seed"),
  VALUE(POP_SIZE, size_t, 512, "Population size."),
  VALUE(MAX_GENS, size_t, 50000, "Maximum number of generations."),
  VALUE(MAX_EVALS, size_t, 2560000000, "Maximum number of trait evaluations."),
  VALUE(STOP_MODE, std::string, "generations", "How do we know when to stop? Options: generations\nevaluations"),

  GROUP(SELECTION, "Selection settings"),
  VALUE(SELECTION, std::string, "tournament", "Selection scheme to use"),
  VALUE(TOURNAMENT_SIZE, size_t, 4, "Tournament size for selection schemes that use tournaments"),

  GROUP(PSB_PROBLEM, "Problem-related settings"),
  VALUE(PROBLEM, std::string, "small-or-large", "Problem to solve"),
  VALUE(TESTING_SET_PATH, std::string, "testing.json", "Path to testing set (json)"),
  VALUE(TRAINING_SET_PATH, std::string, "training.json", "Path to training set (json)"),

  GROUP(POP_INIT, "Population initialization settings"),
  VALUE(POP_INIT_MODE, std::string, "random", "How should we initialize the population? Options: random, load"),
  VALUE(ANCESTOR_FILE_PATH, std::string, "ancestor.json", "Path to ancestor file"),

  GROUP(EVALUATION, "How are organisms evaluated?"),
  VALUE(EVAL_MODE, std::string, "full", "Evaluation mode. Options:\nfull\ncohort\ndown-sample"),
  VALUE(EVAL_FIT_EST_MODE, std::string, "none", "Fitness function estimation method. Options:\nnone\nancestor\nrelative"),
  VALUE(EVAL_MAX_PHYLO_SEARCH_DEPTH, int, -1, "Maximum phylogeny search depth when estimating fitness."),
  VALUE(EVAL_CPU_CYCLES_PER_TEST, size_t, 128, "Maximum number of CPU cycles programs are run for single test case."),
  VALUE(EVAL_ADJ_EST, bool, false, "Penalize estimates by distance?"),
  VALUE(NUM_COHORTS, size_t, 2, "How many cohorts should we divide the tests and organisms into?"),
  VALUE(TEST_DOWNSAMPLE_RATE, double, 0.5, "Proportion of training cases to down-sample each generation"),

  GROUP(SGP_CPU, "SignalGP Virtual CPU"),
  VALUE(MAX_ACTIVE_THREAD_CNT, size_t, 8, "Maximum number of active threads that can run simultaneously on a SGP virtual CPU."),
  VALUE(MAX_THREAD_CAPACITY, size_t, 16, "Maximum thread capacity."),

  GROUP(SGP_PROGRAM, "SignalGP Program settings"),
  VALUE(PRG_MIN_FUNC_CNT, size_t, 0, "Minimum number of functions per program."),
  VALUE(PRG_MAX_FUNC_CNT, size_t, 64, "Maximum number of functions per program."),
  VALUE(PRG_MIN_FUNC_INST_CNT, size_t, 0, "Minimum number of instructions per function."),
  VALUE(PRG_MAX_FUNC_INST_CNT, size_t, 128, "Maximum number of instructions per function."),
  VALUE(PRG_INST_MIN_ARG_VAL, int, -4, "Minimum instruction-argment value"),
  VALUE(PRG_INST_MAX_ARG_VAL, int, 4, "Maximum instruction-argument value"),

  GROUP(SGP_MUTATION, "SignalGP mutation settings"),
  VALUE(MUT_RATE_INST_ARG_SUB, double, 0.005, "InstArgSub rate"),
  VALUE(MUT_RATE_INST_SUB, double, 0.005, "InstSub rate"),
  VALUE(MUT_RATE_INST_INS, double, 0.005, "InstIns rate"),
  VALUE(MUT_RATE_INST_DEL, double, 0.005, "InstDel rate"),
  VALUE(MUT_RATE_SEQ_SLIP, double, 0.05, "SeqSlip rate"),
  VALUE(MUT_RATE_FUNC_DUP, double, 0.05, "FuncDup rate"),
  VALUE(MUT_RATE_FUNC_DEL, double, 0.05, "FuncDel rate"),
  VALUE(MUT_RATE_INST_TAG_BF, double, 0.001, "InstArgTagBF rate"),
  VALUE(MUT_RATE_FUNC_TAG_BF, double, 0.001, "FuncTagBF rate"),
  VALUE(MUT_RATE_INST_TAG_SINGLE_BF, double, 0.0, "Per-tag single bit flip rate"),
  VALUE(MUT_RATE_FUNC_TAG_SINGLE_BF, double, 0.0, "Per-tag single bit flip rate"),
  VALUE(MUT_RATE_INST_TAG_SEQ_RAND, double, 0.0, "Per-tag sequence randomization rate"),
  VALUE(MUT_RATE_FUNC_TAG_SEQ_RAND, double, 0.0, "Per-tag sequence randomization rate"),

  GROUP(OUTPUT, "Output settings"),
  VALUE(OUTPUT_DIR, std::string, "./output/", "What directory are we dumping all this data"),
  VALUE(OUTPUT_SUMMARY_DATA_INTERVAL, size_t, 10, "How often should we output summary data?"),
  VALUE(PRINT_INTERVAL, size_t, 1, "How often do we print run status information?"),
  VALUE(SNAPSHOT_INTERVAL, size_t, 100, "How often should we snapshot?"),
  VALUE(RECORD_PHYLO_GENOTYPES, bool, false, "Output phylogeny genotypes when taking a phylogeny snapshot?"),

  GROUP(MUTATION_ANALYSIS, "Mutation analysis settings"),
  VALUE(MUTATION_ANALYSIS_MODE, bool, false, "Run in mutation analysis mode?"),
  VALUE(NUM_MUTANTS, size_t, 100, "How many mutants to sample?"),
  VALUE(FOCAL_GENOTYPES_FPATH, std::string, "programs.sgp", "Path to genotype file")

)

}